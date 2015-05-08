from __future__ import division
__author__ = 'chris'

import numpy as np
import matplotlib.pyplot as plt
import logging
from psth_and_rasters import get_rasters, make_spike_array
from data_classes import Recording
from odor_analysis import get_odor_name
from utils.h5_decorator import h5decorator


@h5decorator
def find_laser_events(h5, amplitude='', onset_delay=None, *args, **kwargs):
    """

    :param h5:
    :param amplitude: Optional string. Default ''. If specified, include only trials with amplitude matching this value.
    :param onset_delay: None or Int. Default None. If None, include all onset delays; if numeric, include trials with
    matching onset delay.
    :param args:
    :param kwargs:
    :return:
    """
    assert isinstance(h5, Recording)
    laser_events = h5.root.events.laser.events.read()
    laser_params = h5.root.events.laser.params_table.read()
    matches = np.ones(laser_events.shape, dtype=np.bool)
    trial_params = h5.root.events.trials.params_table.read()
    trial_events = h5.root.events.trials.events.read()

    #  filter based on parameter input
    #  TODO: make these work diffently for multiple channels. (ie by taking a tuple for the parameter in)
    if amplitude:
        matches *= laser_params['amplitude_1'] == amplitude
    if onset_delay:
        foo = laser_params['pulseOnsetDelay_1'] == onset_delay
        onset_mask = np.hstack([foo[:, np.newaxis], foo[:, np.newaxis]])
        matches *= onset_mask

    match_idxes = np.where(matches)[0]  # CANNOT USE LOGICAL INDEXING OR THE ARRAY IS FLATTENED!

    return laser_events[match_idxes]


@h5decorator
def find_laser_trains(h5, min_dist_ms=20, onset_delay=None, n_pulses=-1, include_odor_trials=False, *args, **kwargs):
    """
    Finds the start and ends of laser trains. Returns an 2d array where the first column is made of start times and
    2nd column is end times. For each row, the start corresponds to the start of the first pulse in the train, while the
    end corresponds to the offset of the final pulse in the train.

    This uses the find_laser_events method to filter pulses based on amplitude and onset delay if the parameters are
    passed. THIS DOES NOT HANDLE HETEROGENEOUS TRAINS WELL (ie first pulse has some characteristic while other pulses
    are different).


    :param h5:
    :param min_dist_ms: Pulses spaced more than this are considered in separate trains.
    :param include_odor_trials: Bool or string. If false, include no odor trials. If string, include odor trials with
    odors closely matching the string (get_odor_name).
    :param args:
    :param kwargs:
    :return:
    """
    assert isinstance(h5, Recording)
    laser_events = find_laser_events(h5, onset_delay=onset_delay, *args, **kwargs)  # this does not handle the case where laser
    fs = h5.root.events._v_attrs['sample_rate_Hz']
    min_dist = np.int(fs // 1000 * min_dist_ms)
    diffs = np.diff(laser_events[:, 0]) #  pulse starts.
    start_idxes = []  # these indeces correspond to laser event times.
    end_idxes = []  # the pulse after these come at least min_dist away in time.
    num_pulses = []  # the number of pulses that occur between start and end idxes.
    n_pul = 0  # running count
    for i in xrange(0, len(diffs)):
        if diffs[i] < min_dist and diffs[i-1] >= min_dist:
            start_idxes.append(i)
            n_pul += 1
        elif diffs[i] < min_dist and diffs[i-1] < min_dist:
            n_pul += 1

        elif diffs[i] >= min_dist and diffs[i-1] < min_dist:
            n_pul += 1
            end_idxes.append(i)
            num_pulses.append(n_pul)
            n_pul = 0
        elif diffs[i] >= min_dist and diffs[i-1] >= min_dist:
            start_idxes.append(i)
            end_idxes.append(i)
            num_pulses.append(1)
            n_pul = 0
        else:
            raise Exception("This shouldn't happen!!")

    if len(start_idxes) > len(end_idxes):  # if the last pulse was a "start", then it will have no end.
        end_idxes.append(start_idxes[-1])

    if n_pulses > 0:
        num_pulses = np.asarray(num_pulses)
        n_pulses_mask = num_pulses == n_pulses
    else:
        n_pulses_mask = np.ones(len(num_pulses), dtype=np.bool)

    train_events = np.zeros((np.sum(n_pulses_mask), 2), dtype=laser_events.dtype)

    start_idxes = np.asarray(start_idxes)
    end_idxes = np.asarray(end_idxes)

    start_idxes = start_idxes[n_pulses_mask]
    end_idxes = end_idxes[n_pulses_mask]

    train_events[:, 0] = laser_events[start_idxes, 0]
    train_events[:, 1] = laser_events[end_idxes, 1]

    # filter if we don't want to include odor trials or want to filter by odor name.
    if not include_odor_trials or isinstance(include_odor_trials, str):
        new_mask = np.zeros(len(train_events), dtype=np.bool)
        trial_params = h5.root.events.trials.params_table.read()
        trial_evs = h5.root.events.trials.events.read()
        trial_starts = trial_evs[:, 0]
        trial_ends = trial_evs[:, 1]
        if isinstance(include_odor_trials, str):
            odor = get_odor_name(h5, include_odor_trials)

        for i in xrange(len(train_events)):
            train_ev = train_events[i, :]
            train_start = train_ev[0]
            trial_matches = (trial_starts < train_start) * (trial_ends > train_start)
            if not np.sum(trial_matches) == 1:
                new_mask[i] = False
                continue
            t_params = trial_params[trial_matches]
            if isinstance(include_odor_trials, str) and odor == t_params['odor']:
                new_mask[i] = True
            elif isinstance(include_odor_trials, str) and odor != t_params['odor']:
                new_mask[i] = False
            elif not include_odor_trials and (t_params['odor'][0].lower() != 'blank' and t_params['odor'][0] != ''):
                new_mask[i] = False
            elif not include_odor_trials and (t_params['odor'][0].lower() == 'blank' or t_params['odor'][0] == ''):
                new_mask[i] = True
        train_events = train_events[new_mask, :]

    return train_events


@h5decorator
def get_laser_train_rasters(h5, clu, amplitude='', onset_delay=None, n_pulses=-1, time_window_ms=1000, pre_pad_ms=500,
                            include_odor_trials=False, baseline_nsniffs=5, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param amplitude:
    :param onset_delay:
    :param time_window_ms:
    :param pre_pad_ms:
    :return:
    """

    pre_pad_samples = np.int(pre_pad_ms * h5.root.events._v_attrs['sample_rate_Hz'] / 1000.)
    train_events = find_laser_trains(h5, amplitude=amplitude, onset_delay=onset_delay,
                                     n_pulses=n_pulses, include_odor_trials=include_odor_trials, *args, **kwargs)

    train_starts = train_events[:, 0]

    rasters = get_rasters(h5, clu, train_starts-pre_pad_samples, time_window_ms, convert_to_ms=True)
    rasters -= pre_pad_ms

    #  now find a pre - stimulus sniff to calculate baseline from and get a raster centered around this sniff at the
    #  same latency as the laser stimulus.
    baseline_starts = []
    sniff_events = h5.root.events.sniff.inh_events
    sniff_starts = sniff_events[:, 0]
    for start in train_starts:

        last_inh = sniff_starts[sniff_starts < start][-1]
        latency = start - last_inh
        for i in xrange(baseline_nsniffs):
            two_sniffs_back = sniff_starts[sniff_starts < start][-3-i]  # find the inhalation 2 prior to the triggering one.
            baseline_starts.append(two_sniffs_back+latency)  # add the latency, because we want to center around the same point in the sniff as we did in the laser trial.

    baseline_starts = np.asarray(baseline_starts)
    baseline_rasters = get_rasters(h5, clu, baseline_starts - pre_pad_samples, time_window_ms, convert_to_ms=True)
    baseline_rasters -= pre_pad_ms

    return rasters, baseline_rasters


@h5decorator
def plot_laser_train_rasters(h5, clu, amplitude='', onset_delay=None, n_pulses=-1, time_window_ms=1000, pre_pad_ms=500,
                             include_odor_trials=True, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param amplitude:
    :param onset_delay:
    :param time_window_ms:
    :param pre_pad_ms:
    :param args:
    :param kwargs:
    :return:
    """
    rasters, _ = get_laser_train_rasters(h5, clu, amplitude=amplitude, onset_delay=onset_delay, n_pulses=n_pulses,
                                         time_window_ms=time_window_ms, pre_pad_ms=pre_pad_ms,
                                         include_odor_trials=include_odor_trials, *args, **kwargs)
    n_events = len(rasters)
    rows = np.ones(rasters.shape, dtype=np.int)
    row_template = np.arange(n_events)
    rows *= row_template[:, np.newaxis]
    plt.scatter(rasters, rows, marker='|', *args, **kwargs)
    plt.ylim(-1, n_events)
    plt.xlim(-pre_pad_ms, -pre_pad_ms+time_window_ms)
    return


@h5decorator
def get_laser_train_psth(h5, clu, bin_size, amplitude='', onset_delay=None, n_pulses=-1, time_window_ms=1000,
                         pre_pad_ms=500, include_odor_trials=False, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param amplitude:
    :param onset_delay:
    :param time_window_ms:
    :param pre_pad_ms:
    :param args:
    :param kwargs:
    :return: this is returning both the psth and the baseline psths!!
    """
    rstrs = get_laser_train_rasters(h5, clu, amplitude=amplitude, onset_delay=onset_delay, n_pulses=n_pulses,
                                    time_window_ms=time_window_ms, pre_pad_ms=pre_pad_ms,
                                    include_odor_trials=include_odor_trials, *args, **kwargs)
    time_arrays = []
    spike_arrays = []
    n_trialss = []
    psths = []



    for ras in rstrs:
        spike_array, time_array = make_spike_array(ras, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
        spike_arrays.append(spike_array)
        time_arrays.append(time_array)
        n_trials = spike_array.shape[1]
        n_trialss.append(n_trials)
        psth = np.sum(spike_array, axis=1)
        psth_hz = (1000 / bin_size) * psth / n_trials
        psths.append(psth_hz)

    return psths, time_arrays, n_trialss


@h5decorator
def get_laser_train_psth_by_trial(h5, clu, bin_size, amplitude='', onset_delay=None, n_pulses=-1,
                         time_window_ms=1000, pre_pad_ms=500, include_odor_trials=False, *args, **kwargs):
    """

    :param h5:
    :param clu:
    :param amplitude:
    :param onset_delay:
    :param time_window_ms:
    :param pre_pad_ms:
    :param args:
    :param kwargs:
    :return: this is returning both the psth and the baseline psths!!
    """
    rstrs = get_laser_train_rasters(h5, clu, amplitude=amplitude, onset_delay=onset_delay, n_pulses=n_pulses,
                                    time_window_ms=time_window_ms, pre_pad_ms=pre_pad_ms,
                                    include_odor_trials=include_odor_trials, *args, **kwargs)
    time_arrays = []
    spike_arrays = []
    n_trialss = []
    psths = []



    for ras in rstrs:
        spike_array, time_array = make_spike_array(ras, bin_size, -pre_pad_ms, time_window_ms-pre_pad_ms)
        spike_arrays.append(spike_array)
        time_arrays.append(time_array)


    return spike_arrays, time_arrays


@h5decorator
def plot_laser_train_psth(h5, clu, binsize, amplitude='', onset_delay=None, n_pulses=-1,
                          time_window_ms=1000, pre_pad_ms=500, plot_baseline=True,
                          plt_params='', plt_label='', include_odor_trials=False,
                          *args, **kwargs):


    psths, time_arrays, n_trialss = get_laser_train_psth(h5, clu, binsize, amplitude=amplitude, onset_delay=onset_delay,
                                                         n_pulses=n_pulses, time_window_ms=time_window_ms,
                                                         pre_pad_ms=pre_pad_ms, include_odor_trials=include_odor_trials,
                                                         *args, **kwargs)
    styles = ['b', '--k']

    plt.plot(time_arrays[0], psths[0], plt_params, label=plt_label, *args, **kwargs)

    if plot_baseline:
        plt.plot(time_arrays[1], psths[1], '--k', label=plt_label, *args, **kwargs)

    plt.xlim([-pre_pad_ms, -pre_pad_ms+time_window_ms])
    # plt.ylim([0, plt.ylim()[1]])
    plt.ylabel('Firing rate (Hz)')
    plt.xlabel('Time post laser onset (ms)')
    return