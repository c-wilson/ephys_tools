from __future__ import division
__author__ = 'chris'

import tables as tb
import numpy as np
import matplotlib.pyplot as plt
import logging
from psth_and_rasters import get_rasters, make_spike_array


def find_odor_events(h5, odor_name, odor_concentration=None, exact_odor_match=False):
    """
    Finds all event times in which match a specified odor name and concentration.

    :param h5: File object. (This is 'self' when used within a class).
    :param odor_name:
    :param odor_concentration:
    :param
    :return:
    """
    assert isinstance(h5, tb.File)
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
    return fv_events[matches]

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


def find_sniff_events(h5, fv_event):
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
    return sniff_events[ind]


def get_odor_rasters(h5, clu, odor, odor_conc, time_window_ms=1000, pre_pad_ms=500):
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
    fv_events = find_odor_events(h5, odor, odor_conc)
    starts = []
    for i, ev in enumerate(fv_events):
        sn_event = find_sniff_events(h5, ev)[0]  # the first inhalation after odor onset.
        try:
            starts.append(sn_event[0])
        except IndexError:   # empty array, no sniffs found
            print 'no sniffs found'
            pass
    starts = np.array(starts)
    rasters = get_rasters(h5, clu, starts-pre_pad_samples, time_window_ms=time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    return rasters


def plot_odor_rasters(h5, clu, odor, odor_conc, time_window_ms=1000, pre_pad_ms=500, *args, **kwargs):
    """
    Returns array of

    :param h5:
    :param clu:
    :param odor:
    :param odor_conc:
    :return:
    """

    rasters = get_odor_rasters(h5, clu, odor, odor_conc, time_window_ms, pre_pad_ms)
    n_events = len(rasters)
    rows = np.ones(rasters.shape)
    row_template = np.arange(n_events)
    rows *= row_template[:, np.newaxis]
    plt.scatter(rasters, rows, *args, **kwargs)

    return rasters


def get_odor_psth(h5, clu, odor, odor_conc, bin_size, time_window_ms=1000, pre_pad_ms=0):
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
    rasters = get_odor_rasters(h5, clu, odor, odor_conc, time_window_ms, pre_pad_ms)
    spike_array, time_array = make_spike_array(rasters, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
    psth = np.sum(spike_array, axis=1)
    n_trials = spike_array.shape[1]
    return psth, time_array, n_trials

def get_baseline_psth(h5, clu, bin_size, time_window_ms, pre_pad_ms):
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
    fv_events = h5.root.events.finalvalve.events
    starts = []
    for fv_ev in fv_events:
        fv_start = fv_ev[0]
        last_inh = sniff_starts[sniff_starts < fv_start][-1]
        starts.append(last_inh)
    starts = np.array(starts)
    rasters = get_rasters(h5, clu, starts - pre_pad_samples, time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms
    spike_array, time_array = make_spike_array(rasters, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
    psth = np.sum(spike_array, axis=1)
    n_trials = len(starts)
    return psth, time_array, n_trials


def plot_odor_psth_no_baseline(h5, clu, odor, odor_conc, bin_size, time_window_ms=1000, pre_pad_ms=500, *args, **kwargs):
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
    psth, time, _ = get_odor_psth(h5, clu, odor, odor_conc, bin_size, time_window_ms, pre_pad_ms)
    plt.plot(time, psth, *args, **kwargs)
    return

def plot_odor_psth_w_baseline(h5, clu, odor, odor_conc, bin_size, time_window_ms=1000, pre_pad_ms=500, *args, **kwargs):
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
    psth, time, n_odor_trials = get_odor_psth(h5, clu, odor, odor_conc, bin_size, time_window_ms, pre_pad_ms)
    plt.plot(time, psth, *args, **kwargs)
    base_psth, time, n_base_trials = get_baseline_psth(h5, clu, bin_size, time_window_ms, pre_pad_ms)
    base_psth_norm = base_psth * (n_odor_trials/n_base_trials)
    plt.plot(time, base_psth_norm, '--r')
    plt.ylim([0, n_odor_trials])
    return

def plot_odor_psth_baseline_subtract():
    pass

def null(ksndflksdn):
    """



    :param ksndflksdn:
    :return:
    """
    pass




