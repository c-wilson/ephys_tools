from __future__ import division
__author__ = 'chris'

import tables as tb
import numpy as np
from scipy.signal import firwin, filtfilt, find_peaks_cwt
import matplotlib.pyplot as plt
import numba
import logging
import argparse
from utils.detect_peaks import detect_peaks


def make_sniff_events(h5, overwrite=True, *args, **kwargs):
    """

    :param h5: tables.File object.
    :return:
    """
    assert isinstance(h5, tb.File)
    events_group = h5.root.events
    try:
        sniff_stream = h5.root.streams.sniff
    except tb.NoSuchNodeError:
        print('Warning, no sniff stream found.')
        return

    try:
        sniff_events_group = h5.create_group('/events', 'sniff')
    except tb.NodeError as e:
        if overwrite:
            h5.remove_node('/events/sniff', recursive=True)
            h5.flush()
            sniff_events_group = h5.create_group('/events', 'sniff')
        else:
            raise SniffExistExemption

    fs_stream = sniff_stream._v_attrs['sample_rate_Hz']
    fs_events = events_group._v_attrs['sample_rate_Hz']
    sniff_events_array = find_inhalations_simple(sniff_stream)
    sniff_events_array *= (fs_events/fs_stream)
    sniff_events_array = sniff_events_array.astype(np.int)
    sniff_events = h5.create_carray(sniff_events_group, 'inh_events', obj=sniff_events_array)
    sniff_stats = basic_sniff_stats(sniff_events_array, sniff_stream, h5, fs_stream, fs_events)
    return

def basic_sniff_stats(sniff_events_array, sniff_stream, h5, fs_stream, fs_events):
    """

    :param sniff_events: numpy array with inhalation onsets and offsets.
    :param h5: destination file
    :return: tables.Table
    """

    ev_fs_to_stream = fs_events/fs_stream
    assert ev_fs_to_stream % 1. < 1e-2, 'possible rounding error in make_sniff_events indexing'
    ev_fs_to_stream = np.int(np.round(ev_fs_to_stream))

    table_desc = {'inh_integral': tb.Float64Col(pos=0),
                  'inh_length': tb.Float64Col(pos=1),
                  'exh_length': tb.Float64Col(pos=2),
                  'sniff_length': tb.Float64Col(pos=3)}

    sn_tb = h5.create_table('/events/sniff', 'basic_analysis', table_desc, expectedrows=len(sniff_events_array))
    sniff = -sniff_stream[:, 0]

    for i in xrange(len(sniff_events_array)):
        start = sniff_events_array[i, 0]
        end = sniff_events_array[i, 1]
        row = sn_tb.row
        row['inh_integral'] = calc_inh_integral(sniff, start, end, ev_fs_to_stream)
        row['inh_length'] = end - start
        try:
            row['exh_length'] = sniff_events_array[i+1, 0] - end
            row['sniff_length'] = sniff_events_array[i+1, 0] - start
        except IndexError:
            row['exh_length'] = np.nan
        row.append()
    sn_tb.flush()
    return sn_tb


@numba.autojit('f8(f8[:], i8, i8, i8)') # numba gives ~50-fold performance increase for 200 ms inhalations.
def calc_inh_integral(sniff, start, stop, ev_to_stream_factor):
    """
    calculate the integral of the sniff between two values.
    :param sniff: np.array of sniff stream
    :param start: start sample
    :param stop: stop sample
    :param ev_to_stream_factor: resample factor if start and stop events and stream are not on same time base.
    Typically 16.
    :return:
    """

    integral = 0.
    start = start // ev_to_stream_factor
    stop = stop // ev_to_stream_factor
    for i in xrange(start, stop):
        integral += sniff[i]
    return integral



def find_inhalations_simple(sniff_array, *args, **kwargs):
    """

    :param sniff_array:
    :return:
    """

    assert(isinstance(sniff_array, tb.Array))
    sniff = -sniff_array[:, 0]
    sniff -= np.mean(sniff)
    fs = sniff_array.get_attr('sample_rate_Hz')
    sniff = filt(sniff, fs, cutoff=120, n_taps=41)
    top = np.percentile(sniff, 99)
    bottom = np.percentile(sniff, 1)
    peaks_inh = detect_peaks(sniff,
                             mph=top/10,
                             mpd=(0.1 * fs),  # min peak distance is 100 ms.
                             edge='rising')    
    peaks_exh = detect_peaks(-sniff,
                             mph=bottom/5,
                             mpd=(0.1 * fs),
                             edge='rising')
    sniff_logical = sniff > top/20
    edges = np.convolve([1, -1], sniff_logical, mode='same')

    #quick filter for extremely short sniffs:

    edges_filt = np.zeros(edges.shape, dtype=np.int8)
    min_inh_width_ms = 25
    min_inh_width_samp = np.int(min_inh_width_ms/1000 * fs)
    reject_down = -1
    widths = []
    for i in xrange(edges.size):
        reject = False
        if edges[i] == 1:
            first_pk = (peaks_inh > i)[0]
            first_valley = (peaks_exh > i)[0]
            ex_range = i+3000  # look for an exhalation in the next 3000 samples
            if ex_range > edges.size:  # unless there are less than 3000 samples in the recording remaining.
                ex_range = edges.size
            for ii in xrange(i, ex_range):
                if edges[ii] == -1:
                    next_exh = ii
                    break
            widths.append(next_exh - i)
            if (next_exh - i) <= min_inh_width_samp:
                reject = True
                reject_down = next_exh
                logging.debug('rejection of inh at sample (width): {0}'.format(i))
            elif (first_pk - i) > (first_valley - i):
                reject = True
                reject_down = next_exh
                logging.debug('rejection of inh at sample (peaks): {0}'.format(i))
            else:
                pass

            if not reject:
                edges_filt[i] = 1
            else:
                pass
        elif edges[i] == -1:
            if i != reject_down:
                edges_filt[i] = -1
            else:
                logging.debug('rejection of exh at sample (peaks): {0}'.format(i))

    inh_onsets = np.where(edges_filt == 1)[0]
    inh_offsets = np.where(edges_filt == -1)[0]
    
    # first inhalation should not start on first sample (we missed the actual start, so this isn't really a sniff):
    if inh_onsets[0] == 0:
        inh_onsets = inh_onsets[1:]
        inh_offsets = inh_offsets[1:]
    # first inhalation should be before first exhalation:
    if inh_offsets[0] < inh_onsets[0]:
        inh_offsets = inh_offsets[1:]
    # the arrays need to be the same size, so crop the bigger one.
    while inh_onsets.size > inh_offsets.size:
        inh_onsets = inh_onsets[:-1]
    while inh_onsets.size < inh_offsets.size:
        inh_offsets = inh_offsets[:-1]


    # Bug check, on doesn't come before the next off:
    for i in xrange(inh_onsets.size - 1):
        on = inh_onsets[i]
        n_on = inh_onsets[i+1]
        off = inh_offsets[i]
        assert n_on > off, 'Two consecutive inhalation onsets without an offset were found. Bug.'
        assert on < off, 'Exhalation occurs before inhalation. Bug.'

    #return events array:
    return np.array([inh_onsets, inh_offsets]).T


def filt(x, fs, cutoff=120., n_taps=41):
    """
    Zero delay, linear phase lowpass filter for sniff.

    :param x: sniff array
    :param fs: sampling frequency of x.
    :param cutoff: cut frequency.
    :param n_taps: filter order+1. This should be odd to produce a Type I filter.
    :return:
    """

    a = 1.
    b = firwin(n_taps, cutoff/fs, window='hamming')

    y = filtfilt(b, a, x, axis=0)
    return y


class SniffExistExemption(Exception):
    pass

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description='Run sniff processing.')

    # TODO:
    parser.add_argument('input_filename',
                        help='.prm or .kwik filename')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='only run clustering on first few seconds of edata')
    parser.add_argument('--append', action='store_true', default=False,
                        help='append data to the destination file if kwik file is modified since destination file was created')
    parser.add_argument('--overwrite', action='store_true', default=False,
                       help='overwrite the destination file if it exists')
    args = parser.parse_args()

    with tb.open_file(args.input_filename, 'a') as fi:
        logging.info('starting sniff processing')
        make_sniff_events(fi)
        logging.info('complete.')