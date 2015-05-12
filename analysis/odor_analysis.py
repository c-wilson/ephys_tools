from __future__ import division
__author__ = 'chris'

import tables as tb
import numpy as np
import matplotlib.pyplot as plt
import logging
from psth_and_rasters import get_rasters, make_spike_array
from utils.h5_decorator import h5decorator

@h5decorator
def find_odor_events(h5, odor_name, odor_concentration=None, exact_odor_match=False, include_laser_trials=False,
                     only_laser_trials=False, skip_last_trials=0, *args, **kwargs):
    """
    Finds all event times in which match a specified odor name and concentration.

    :param h5: File object. (This is 'self' when used within a class).
    :param odor_name:
    :param odor_concentration:
    :param
    :return:
    """

    fv_events = h5.root.events.finalvalve.events.read()
    fv_params = h5.root.events.finalvalve.params_table.read()
    odor_matches = fv_params['odor'] == get_odor_name(h5, odor_name, not exact_odor_match)
    if odor_concentration is not None:
        conc_matches = fv_params['odorconc'] == odor_concentration
        logging.debug('for concentration {0}, {1} matches found in FV params.'.format(odor_concentration, conc_matches))
        matches = odor_matches * conc_matches
    else:
        matches = odor_matches
    logging.debug('Total matches found for odor:conc pair: {0}.'.format(len(odor_matches)))

    odor_matches = fv_events[matches]

    laser_trials = np.zeros(len(odor_matches), dtype=np.bool)
    if not include_laser_trials or only_laser_trials:
        # TODO: this is a stupid hack for pre- 24 March recordings where the sniff could saturate the laser channel.
        try:
            laser_params = h5.root.events.trials.params_table[:]['amplitude_1']
            laser_trial_events = h5.root.events.trials.events[:]
            for i in xrange(len(odor_matches)):
                trial_log = (odor_matches[i, 0] > laser_trial_events[:, 0]) * \
                            (odor_matches[i, 0] < laser_trial_events[:, 1])
                assert np.sum(trial_log) == 1
                trial_idx = np.where(trial_log)[0]
                if laser_params[trial_idx]:
                    laser_trials[i] = True
            if only_laser_trials:
                odor_matches = odor_matches[laser_trials]
            else:
                odor_matches = odor_matches[np.invert(laser_trials)]
        except ValueError as e:
            #  if we're missing an amplitude for the laser, than we don't care about laser trials.
            #  if we do have this field in our h5, then the ValueError is something interesting that we want to know.
            if 'amplitude_1' in h5.root.events.trials.params_table.dtype.names:
                raise e
            else:
                pass
    end = len(odor_matches) - skip_last_trials
    return odor_matches[:end]


@h5decorator
def get_odor_name(h5, odor, close_match=True):
    """
    Find close matches to given odor name within file parameters table, correcting for case and partial matches.

    Examples:
    odor strings '(+)-alpha-pinene' and '(+)-limonene' are in the h5 file.
    1) Case: '(+)-alpha-Pinene' will be returned as '(+)-alpha-pinene'
    2) Partial matches: 'pinene' will be returned as '(+)-alpha-pinene'

    If the partial matches are found for two unique strings in the h5 file, an exemption will be raised (ie the string
    'ene' is found in both '(+)-alpha-pinene' and '(+)-limonene').

    :param h5:
    :param odor: string to be checked and corrected if possible.
    :return:
    """

    if not close_match:
        return odor

    odor_params = h5.root.events.finalvalve.params_table[:]['odor']
    if odor not in odor_params:
        matches = []
        odor = odor.lower()
        for o in odor_params:
            if odor in o.lower() or o.lower() in odor:
                matches.append(o)
            # check that all of the matches are the same, so we know that we don't have multiple matches (ie a single
            # letter may not be specific enough, as it can be in multiple different strings).
        if matches:
            if matches[1:] == matches[:-1]:
                logging.info('Input odor "{0}" not found, using close match "{1}".'.format(odor, matches[0]))
                odor = matches[0]
            else:
                multiple_match_set = set(matches)
                raise ValueError('Multiple matches found for input odor "{0}". Matches found in finalvalve paramaters:{1}'.format(odor, multiple_match_set))
        else:
            raise ValueError('odor "{0}" not found in finalvalve parameters.'.format(odor))
    return odor


@h5decorator
def get_unique_odors(h5, concentration=None):
    """
    Returns all unique odor strings found in finalvalve events parameter table. Can be masked by concentration (ie to
    find all of the odors that were presented at a specific concentration or concentrations).


    :param h5:
    :param concentration: can be a list or scalar used to mask trial parameters. Default is None.
    :return:
    """

    params = h5.root.events.finalvalve.params_table[:]
    if concentration:
        mask = np.in1d(params['odorconc'], concentration, assume_unique=True)
        # in1d returns boolean matrix where True indicates that a value in array a is also in b.
        params = params[mask]
    return np.unique(params['odor'])


@h5decorator
def get_unique_concentrations(h5, odor=None, close_odor_match=True):
    """
    Returns all unique odor concentrations found in finalvalve events parameter table. Can be masked by odor (ie to find
    all of the concentrations that were presented for a specific odor or a list of odors).


    :param h5:
    :param odor: string or list of strings used to mask trial parameters. Default is none.
    :param close_odor_match: should we try to find close odor matches, or should we only use exact matches?
    :return:
    """
    params = h5.root.events.finalvalve.params_table[:]
    if odor:
        try:
            odor = get_odor_name(h5, odor)
            mask = np.in1d(params['odor'], odor, assume_unique=True)
            params = params[mask]
        except ValueError as e:
            print e.message
            return []
    return np.unique(params['odorconc'])


@h5decorator
def find_sniff_events(h5, fv_event, *args, **kwargs):
    """
    Finds the first inhalation event within the final valve event.

    :param h5: File object, assumed when used within class structure.
    :param fv_event: event array ([on, off])
    :return:
    """
    sniff_grp = h5.root.events.sniff
    sniff_events = sniff_grp.inh_events[:]
    sniff_starts = sniff_events[:, 0]
    ind = (sniff_starts > fv_event[0])*(sniff_starts < fv_event[1])
    if ind.any():
        return sniff_events[ind]
    else:
        return None


@h5decorator
def get_odor_rasters(h5, clu, odor, odor_conc,
                     time_window_ms=1000,
                     pre_pad_ms=500,
                     include_laser_trials=False,
                     only_laser_trials=False,
                     skip_last_trials=0,
                     *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """
    pre_pad_samples = np.int(pre_pad_ms * h5.root.events._v_attrs['sample_rate_Hz'] / 1000.)
    fv_events = find_odor_events(h5, odor, odor_conc,
                                 include_laser_trials=include_laser_trials,
                                 only_laser_trials=only_laser_trials,
                                 skip_last_trials=skip_last_trials,
                                 *args, **kwargs)
    starts = []
    for i, ev in enumerate(fv_events):
        sn_events = find_sniff_events(h5, ev)
        try:
            sn_event = sn_events[0]  # the first inhalation after odor onset.
            starts.append(sn_event[0])
        except TypeError:   # empty array, no sniffs found
            print 'no sniffs found'
            pass
    starts = np.array(starts)
    rasters = get_rasters(h5, clu, starts-pre_pad_samples, time_window_ms=time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    return rasters


@h5decorator
def get_odor_rasters_sniff(h5, clu, odor, odor_conc, time_window_ms=1000, pre_pad_ms=500, include_laser_trials=False,
                           only_laser_trials=False, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """

    fs = h5.root.events._v_attrs['sample_rate_Hz']
    pre_pad_samples = np.int(pre_pad_ms * fs / 1000.)
    fv_events = find_odor_events(h5, odor, odor_conc, include_laser_trials=include_laser_trials,
                                 only_laser_trials=only_laser_trials)
    inh_starts = list()
    sniff_events_all_trials = list()
    for i, ev in enumerate(fv_events):
        sn_events = find_sniff_events(h5, ev, *args, **kwargs)
        try:
            first_sn_event = sn_events[0]  # the first sniff after odor onset.
            first_inh = first_sn_event[0]  # this is the first inhalation after odor onset. THIS IS TIME 0!!!
            inh_starts.append(first_inh)

            sn_events_in_window = find_sniff_events(h5, [first_inh - pre_pad_samples,
                                                         first_inh - pre_pad_samples + (time_window_ms*fs)/1000.])
            sn_events_in_window = sn_events_in_window - first_inh  # must make relative to first inh = 0.
            sn_events_in_window_ms = sn_events_in_window / fs * 1000
            sniff_events_all_trials.append(sn_events_in_window_ms)
        except TypeError:   # empty array, no sniffs found
            print 'no sniffs found'
            pass

    inh_starts = np.array(inh_starts)
    rasters = get_rasters(h5, clu, inh_starts-pre_pad_samples, time_window_ms=time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    assert len(rasters) == len(sniff_events_all_trials)

    return rasters, sniff_events_all_trials


@h5decorator
def plot_odor_rasters(h5, clu, odor, odor_conc,
                      time_window_ms=1000,
                      pre_pad_ms=500,
                      include_laser_trials=False,
                      only_laser_trials=False,
                      axis=None,
                      skip_last_trials=0,
                      quick_plot=True,
                      *args, **kwargs):
    """
    Returns array of

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :return:
    """
    rasters = get_odor_rasters(h5, clu, odor, odor_conc, time_window_ms, pre_pad_ms,
                               include_laser_trials=include_laser_trials, 
                               only_laser_trials=only_laser_trials,
                               skip_last_trials=skip_last_trials)
    n_events = len(rasters)
    rows = np.ones(rasters.shape)
    row_template = np.arange(n_events) +1
    rows *= row_template[:, np.newaxis]
    if axis is None:
        axis = plt.axes()

    if quick_plot:
        axis.scatter(rasters, rows, marker='|', *args, **kwargs)
        axis.set_ylim(0, n_events+1)
    else:
        for tr, tr_ras in enumerate(rasters):
            axis.vlines(tr_ras, tr+.5, tr+1.5, *args, **kwargs)
            axis.spines['top'].set_visible(False)
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_visible(False)
            axis.xaxis.set_ticks_position('bottom')
            axis.set_yticks([])


    return rasters


@h5decorator
def plot_odor_rasters_sniff(h5, clu, odor, odor_conc, 
                            time_window_ms=1000, 
                            pre_pad_ms=500, 
                            include_laser_trials=False,
                            only_laser_trials=False, 
                            axis=None, 
                            skip_last_trials=0,
                            *args, **kwargs):
    """
    Returns psth for all trials in which odor and concentration are met. Raster time 0 is the start of the first
    inhalation following odor onset.

    PSTH is not normalized or smoothed (ie it is simply a sum of the total number of spikes that for each bin across
    all odor events in the h5.

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """
    rasters, sniff_by_trial = get_odor_rasters_sniff(h5, clu, odor, odor_conc, time_window_ms, pre_pad_ms,
                                                     include_laser_trials=include_laser_trials,
                                                     only_laser_trials=only_laser_trials, 
                                                     skip_last_trials=skip_last_trials,
                                                     *args, **kwargs)
    n_events = len(rasters)
    rows = np.ones(rasters.shape)
    row_template = np.arange(n_events)
    rows *= row_template[:, np.newaxis]
    if axis is None:
        axis = plt.axes()
    axis.scatter(rasters, rows, marker='|', *args, **kwargs)

    for i, trial in enumerate(sniff_by_trial):
        for sniff in trial:
            if sniff[1] > time_window_ms-pre_pad_ms:
                sniff[1] = time_window_ms - pre_pad_ms
            axis.plot(sniff, [i]*2, 'g', linewidth=5, alpha=.2)
    axis.set_ylim(-1, n_events)
    return rasters, sniff_by_trial


@h5decorator
def get_odor_psth(h5, clu, odor, odor_conc, bin_size,
                  time_window_ms=1000,
                  pre_pad_ms=0,
                  skip_last_trials=0,
                  *args, **kwargs):
    """
    Returns psth for all trials in which odor and concentration are met. Raster time 0 is the start of the first
    inhalation following odor onset.

    PSTH is not normalized or smoothed (ie it is simply a sum of the total number of spikes that for each bin across
    all odor events in the h5.

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """

    rasters = get_odor_rasters(h5, clu, odor, odor_conc,
                               time_window_ms=time_window_ms,
                               pre_pad_ms=pre_pad_ms,
                               skip_last_trials=skip_last_trials,
                               *args, **kwargs)
    spike_array, time_array = make_spike_array(rasters, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
    psth = np.sum(spike_array, axis=1)
    n_trials = spike_array.shape[1]
    psth_hz = (1000 / bin_size) * psth / n_trials
    return psth_hz, time_array, n_trials


@h5decorator
def get_pre_odor_baseline_psth(h5, clu, odor, bin_size, n_sniffs=5,
                               include_laser_trials=False,
                               time_window_ms=1000,
                               pre_pad_ms=0,
                               skip_last_trials=0):
    """
    Gets the first inhalations before every finalvalve opening and averages them all to get a baseline psth.

    :param h5:
    :param clu:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """
    pre_pad_samples = np.int(pre_pad_ms * h5.root.events._v_attrs['sample_rate_Hz'] / 1000.)
    sniff_events = h5.root.events.sniff.inh_events
    sniff_starts = sniff_events[:, 0]
    fv_events = find_odor_events(h5, odor,
                                 include_laser_trials=include_laser_trials,
                                 skip_last_trials=skip_last_trials)
    starts = []
    for fv_ev in fv_events:
        fv_start = fv_ev[0]
        for i_sn in xrange(1, n_sniffs + 1):
            last_inh = sniff_starts[sniff_starts < fv_start][-i_sn]
            starts.append(last_inh)
    starts = np.array(starts)
    rasters = get_rasters(h5, clu, starts - pre_pad_samples, time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    spike_array, time_array = make_spike_array(rasters, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
    psth = np.sum(spike_array, axis=1)
    n_trials = len(starts)
    psth_hz = (1000 / bin_size) * psth / len(starts)
    return psth_hz, time_array, n_trials


@h5decorator
def get_pre_odor_baseline_psths(h5, clu, odor, bin_size, n_sniffs=5,
                                include_laser_trials=False,
                                time_window_ms=1000,
                                pre_pad_ms=0,
                                skip_last_trials=0):
    """
    Gets the first inhalations before every finalvalve opening and averages them all to get a baseline psth.

    :param h5:
    :param clu:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """
    pre_pad_samples = np.int(pre_pad_ms * h5.root.events._v_attrs['sample_rate_Hz'] / 1000.)
    sniff_events = h5.root.events.sniff.inh_events
    sniff_starts = sniff_events[:, 0]
    fv_events = find_odor_events(h5, odor,
                                 include_laser_trials=include_laser_trials,
                                 skip_last_trials=skip_last_trials)
    starts = []
    for fv_ev in fv_events:
        fv_start = fv_ev[0]
        for i_sn in xrange(1, n_sniffs + 1):
            last_inh = sniff_starts[sniff_starts < fv_start][-i_sn]
            starts.append(last_inh)
    starts = np.array(starts)
    rasters = get_rasters(h5, clu, starts - pre_pad_samples, time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    spike_array, time_array = make_spike_array(rasters, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
    # psths = np.sum(spike_array, axis=1)
    n_trials = len(starts)
    psths_hz = (1000 / bin_size) * spike_array
    return psths_hz, time_array, n_trials


@h5decorator
def plot_odor_psth_no_baseline(h5, clu, odor, odor_conc, bin_size, include_laser_trials=False, only_laser_trials=False,
                              time_window_ms=1000, pre_pad_ms=500, axis=None, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :param args:
    :param kwargs:
    :return:
    """
    psth, time, _ = get_odor_psth(h5, clu, odor, odor_conc, bin_size,
                                  include_laser_trials=include_laser_trials,
                                  only_laser_trials=only_laser_trials,
                                  time_window_ms=time_window_ms,
                                  pre_pad_ms=pre_pad_ms)
    if axis is None:
        axis = plt.axes()

    axis.plot(time, psth, *args, **kwargs)
    return


@h5decorator
def plot_odor_psth_w_baseline(h5, clu, odor, odor_conc, bin_size,
                              include_laser_trials=False,
                              only_laser_trials=False,
                              time_window_ms=1000,
                              pre_pad_ms=500,
                              axis=None,
                              skip_last_trials=0,
                              *args, **kwargs):
    """



    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :param bin_size:
    :param time_window_ms:
    :param pre_pad_ms:
    :param args:
    :param kwargs:
    :return:
    """

    psth, time, n_odor_trials = get_odor_psth(h5, clu, odor, odor_conc, bin_size,
                                              include_laser_trials=include_laser_trials,
                                              only_laser_trials=only_laser_trials,
                                              time_window_ms=time_window_ms,
                                              pre_pad_ms=pre_pad_ms,
                                              skip_last_trials=skip_last_trials)

    base_psth, time, n_base_trials = get_pre_odor_baseline_psth(h5, clu, odor, bin_size,
                                                                include_laser_trials=include_laser_trials,
                                                                time_window_ms=time_window_ms,
                                                                pre_pad_ms=pre_pad_ms,
                                                                skip_last_trials=skip_last_trials)

    # base_psth_norm = base_psth * (n_odor_trials/n_base_trials)
    if axis is None:
        axis = plt.axes()
    axis.plot(time, psth, *args, **kwargs)
    axis.plot(time, base_psth, '--k')
    axis.set_xlim([-pre_pad_ms, -pre_pad_ms+time_window_ms])
    axis.set_ylim([0, axis.get_ylim()[1]])
    axis.set_ylabel('Firing rate (Hz)')
    axis.set_xlabel('t (ms)')
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.xaxis.set_ticks_position('bottom')
    axis.yaxis.set_ticks_position('left')

    axis.plot([0, 0], axis.get_ylim(), '-k')

    return




