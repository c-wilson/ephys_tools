from __future__ import division
__author__ = 'chris'


import numpy as np
import tables as tb
from scipy import signal
import logging


"""
parses and adds stream edata to overall h5 object.

to add a stream to the system, add to the stream_defs dictionary:
stream_defs = {'YOUR_STREAM_ARRAY_NAME': decimation factor)
"""

stream_defs = {'sniff': 16}


def streamer(raw_kwd_fn, destination):
    logging.info('Adding streams')
    assert isinstance(destination, tb.File)
    with tb.open_file(raw_kwd_fn, 'r') as raw_kwd:
        for s_name, downsample_factor in stream_defs.iteritems():
            logging.info('Adding stream: {0}'.format(s_name))
            add_stream(s_name, raw_kwd, destination,  downsample_factor)
        add_lfp_stream(raw_kwd, destination, downsample_factor)
    logging.info('Streamer complete.')


def add_stream(stream_name, raw_kwd, dest_file, downsample_factor=16, *args, **kwargs):
    """
    Add decimated analog stream data to destination file.

    :param raw_kwd: origin file.
    :param dest_file: destination hdf5 file
    :param stream_name: Name of the stream to be processed (should be the name of the HDF5 array)
    :param downsample_factor: decimate by this factor (ie take every nth sample)
    :param args:
    :param kwargs:
    :return:
    """

    assert isinstance(dest_file, tb.File)
    assert isinstance(raw_kwd, tb.File)
    assert isinstance(downsample_factor, int)

    n_recs = raw_kwd.root.recordings._v_nchildren
    dec_stream = np.array([], dtype=np.int16)
    dec_stream = dec_stream[:, np.newaxis]
    for i in xrange(n_recs):
        rec = raw_kwd.get_node('/recordings/{0:d}'.format(i))
        try:
            s_pointer = raw_kwd.get_node('/recordings/{0:d}/{1:s}'.format(i, stream_name))
            data_len = rec.data.shape[0]  # in case data is truncated relative to the stream.
            s_mem = s_pointer[:data_len]
            s_dec = decimate(s_mem)
            dec_stream = np.concatenate((dec_stream, s_dec))

        except tb.NoSuchNodeError:
            if i == 0:
                logging.warning('No array exists for stream "{0}" in first recording. This is probably because '
                                'this stream was not recorded, so stream will be skipped.'.format(stream_name, i))
            else:
                raise ValueError('No array exists for stream "{0}" in rec {1}.'.format(stream_name, i))
            return
    try:
        fs = s_pointer.get_attr('sampling_rate_Hz')
    except AttributeError:
        logging.warning('No sampling rate specified in {0} stream, using default of 20833.'.format(stream_name))
        fs = 20833.
    st_arr = dest_file.create_carray('/streams',
                                     stream_name,
                                     obj=dec_stream,
                                     createparents=True)
    st_arr.set_attr('sample_rate_Hz', fs/downsample_factor)
    dest_file.flush()
    logging.info('complete.')
    return

def add_lfp_stream(raw_kwd, dest_file, downsample_factor=16, *args, **kwargs):
    """
    Adds decimated LFP stream to destination h5 file.

    :param raw_kwd:
    :param dest_file:
    :param stream_name:
    :param downsample_factor:
    :param args:
    :param kwargs:
    :return:
    """
    logging.info('Adding LFP.')
    assert isinstance(dest_file, tb.File)
    assert isinstance(raw_kwd, tb.File)
    assert isinstance(downsample_factor, int)

    l = 0
    records = []
    nrecs = raw_kwd.root.recordings._v_nchildren
    for i in xrange(nrecs):
        record = raw_kwd.get_node('/recordings/{0:d}/data'.format(i))  # using PL filtered data.
        l += record.shape[0]

        records.append(record)
    n_ch = record.shape[1]
    try:
        o_fs = record.get_attr('sample_rate_Hz')
    except AttributeError:
        o_fs = 20833.
    n_fs = o_fs / downsample_factor
    lfp = dest_file.create_earray('/streams',
                                  'LFP',
                                  atom=tb.Int16Atom(),
                                  shape=(0, n_ch),  # extending in the n_samples direction.
                                  title="local field potential ({0:0.1f} Hz)".format(n_fs),
                                  createparents=True,
                                  filters=tb.Filters(complevel=6))
    lfp.set_attr('sample_rate_Hz', n_fs)
    for i, record in enumerate(records):
        #calculate size of decimated signal:
        s = (record.shape[0] + downsample_factor - 1 - 25) // downsample_factor
        # this is hardcoded with first sample of 25, which corresponds to a decimation n = 50.
        temp_sig = np.zeros((s, n_ch), dtype=np.int16)
        for ii, sig in enumerate(record):
            logging.info('processing LFP in rec {0}, ch {1}'.format(i,ii))
            temp_sig[:, ii] = decimate(sig, downsample_factor)
        lfp.append(temp_sig)
    lfp.flush()
    dest_file.flush()
    logging.info('complete.')
    return



def decimate(x, decimation_factor=16):
    """
    This reduces the sampling rate (frequency) of a signal. It uses a delay-compensated linear phase filter for
    antialiasing. The implementation is the same as that used by SpikeDetect2.

    :param x: array to decimate.
    :param decimation_factor: divide sampling rate by this number (int)
    :return:
    :type decimation_factor: int
    """
    assert isinstance(decimation_factor, int)
    q = decimation_factor
    n = 50  # filter order.

    # using firwin because it gives a linear phase filter (delay for all phases is equal: no "delay distortion")
    b = signal.firwin(n+1, 1. / q, window='hamming')
    a = 1.

    y = signal.lfilter(b, a, x)

    sl = [slice(None)] * y.ndim
    sl[0] = slice(n // 2, None, q)
    # Phase delay for FIR filter is equal to  .5 * (num_taps-1)/fs.
    # Since fs = 1 here, this means that the slicing starts at exactly the offset value, which is the desired behavior.

    return y[sl]



if __name__ == '__main__':
    pass