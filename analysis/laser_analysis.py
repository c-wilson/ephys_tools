from __future__ import division
__author__ = 'chris'

import numpy as np
import matplotlib.pyplot as plt
import logging
from psth_and_rasters import get_rasters, make_spike_array
from data_classes import Recording

def find_laser_events(h5, amplitude='', onset_delay=None, *args, **kwargs):
    """

    :param h5:
    :param amplitude:
    :param onset_delay:
    :param args:
    :param kwargs:
    :return:
    """
    assert isinstance(h5, Recording)
    laser_events = h5.root.events.laser.events.read()
    laser_params = h5.root.events.laser.params_table.read()
    matches = np.ones(laser_events.shape, dtype=np.bool)

    #  filter based on parameter input
    #  TODO: make these work diffently for multiple channels. (ie by taking a tuple for the parameter in)
    if amplitude:
        matches *= laser_params['amplitude_1'] == amplitude
    if onset_delay:
        assert isinstance(onset_delay, int)
        matches *= laser_params['pulseOnSetDelay_1'] == onset_delay

    match_idxes = np.where(matches)[0]  # CANNOT USE LOGICAL INDEXING OR THE ARRAY IS FLATTENED!

    return laser_events[match_idxes]


def find_laser_trains(h5, min_dist_ms=20, *args, **kwargs):
    """
    Finds the start and ends of laser trains. Returns an 2d array where the first column is made of start times and
    2nd column is end times. For each row, the start corresponds to the start of the first pulse in the train, while the
    end corresponds to the offset of the final pulse in the train.

    This uses the find_laser_events method to filter pulses based on amplitude and onset delay if the parameters are
    passed. THIS DOES NOT HANDLE HETEROGENEOUS TRAINS WELL (ie first pulse has some characteristic while other pulses
    are different).


    :param h5:
    :param min_dist_ms: Pulses spaced more than this are considered in separate trains.
    :param args:
    :param kwargs:
    :return:
    """

    assert isinstance(h5, Recording)
    laser_events = find_laser_events(h5, *args, **kwargs)  # this does not handle the case where laser
    fs = h5.root.events._v_attrs['sample_rate_Hz']
    min_dist = np.int(fs // 1000 * min_dist_ms)
    diffs = np.diff(laser_events[:, 0]) #  pulse starts.
    start_idxes = [0]  # the first idx is, by definition the start of a 'train'
    end_idxes = []  # the pulse after these come at least min_dist away in time.
    num_pulses = []  # the number of pulses that occur between start and end idxes.
    n_pulses = 0  # running count
    for i in xrange(1, len(diffs)):
        if diffs[i] < min_dist and diffs[i-1] > min_dist:
            start_idxes.append(i)
            n_pulses += 1
        elif diffs[i] < min_dist and diffs[i-1] < min_dist:
            n_pulses += 1

        elif diffs[i] > min_dist and diffs[i-1] < min_dist:
            n_pulses += 1
            end_idxes.append(i)
            num_pulses.append(n_pulses)
            n_pulses = 0
        elif diffs[i] > min_dist and diffs[i-1] > min_dist:
            start_idxes.append(i)
            end_idxes.append(i)
            num_pulses.append(1)
            n_pulses = 0
        else:
            raise Exception("This shouldn't happen!!")

    if len(start_idxes) > len(end_idxes):  # if the last pulse was a "start", then it will have no end.
        end_idxes.append(start_idxes[-1])

    train_events = np.zeros((len(start_idxes), 2), dtype=laser_events.dtype)
    train_events[:, 0] = laser_events[start_idxes, 0]
    train_events[:, 1] = laser_events[end_idxes, 1]

    return train_events


def get_laser_train_rasters(h5, clu, amplitude='', onset_delay=None, time_window_ms=1000, pre_pad_ms=500):
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
    train_events = find_laser_trains(h5, amplitude=amplitude, onset_delay=onset_delay)

    train_starts = train_events[:,0]

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
        two_sniffs_back = sniff_starts[sniff_starts < start][-3]  # find the inhalation 2 prior to the triggering one.
        baseline_starts.append(two_sniffs_back+latency)  # add the latency, because we want to center around the same point in the sniff as we did in the laser trial.

    baseline_starts = np.asarray(baseline_starts)
    baseline_rasters = get_rasters(h5, clu, baseline_starts - pre_pad_samples, time_window_ms, convert_to_ms=True)
    baseline_rasters -= pre_pad_ms

    return rasters, baseline_rasters


def plot_laser_train_rasters(h5, clu, amplitude='', onset_delay=None, time_window_ms=1000, pre_pad_ms=500,
                             *args, **kwargs):
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
    rasters, _ = get_laser_train_rasters(h5, clu, amplitude, onset_delay, time_window_ms, pre_pad_ms)
    n_events = len(rasters)
    rows = np.ones(rasters.shape, dtype=np.int)
    row_template = np.arange(n_events)
    rows *= row_template[:, np.newaxis]
    plt.scatter(rasters, rows, marker='|', *args, **kwargs)
    plt.plot(rasters)
    plt.ylim(-1, n_events)
    return

def get_laser_train_psth(h5, clu, binsize, amplitude='', onset_delay=None, time_window_ms=1000, pre_pad_ms=500,
                   *args, **kwargs):
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
    rstrs = get_laser_train_rasters(h5, clu, amplitude, onset_delay, time_window_ms, pre_pad_ms)
    time_arrays = []
    spike_arrays = []
    n_trialss = []
    psths = []

    for ras in rstrs:
        spike_array, time_array = make_spike_array(ras, binsize, -pre_pad_ms, time_window_ms-pre_pad_ms)
        spike_arrays.append(spike_array)
        time_arrays.append(time_array)
        psths.append(np.sum(spike_array, axis=1))
        n_trialss.append(spike_array.shape[1])


    return psths, time_arrays, n_trialss

def plot_laser_train_psth(h5, clu, binsize, amplitude='', onset_delay=None, time_window_ms=1000, pre_pad_ms=500,
                          *args, **kwargs):


    psths, time_arrays, n_trialss = get_laser_train_psth(h5, clu, binsize, amplitude, onset_delay,
                                                         time_window_ms, pre_pad_ms)
    styles = ['b', '--r']
    for psth, time, style in zip(psths, time_arrays, styles):
        plt.plot(time, psth, style, *args, **kwargs)
    plt.xlim([-pre_pad_ms, -pre_pad_ms+time_window_ms])
    return