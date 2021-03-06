from __future__ import division
__author__ = 'chris'

import numpy as np
import numba
from utils.h5_decorator import h5decorator


@h5decorator
def get_rasters(h5, clu, start_sample_array, time_window_ms, convert_to_ms=False):
    """
    Return spikes within windows from starts to starts+time_window_ms
    Spike times are returned in samples RELATIVE TO START!!!

    :param clu: Cluster Array or string corresponding to the location of the cluster array.
    :param starts: List or array of where to extract rasters from, in SAMPLES, just like the cluster arrays.
    :param time_window_ms: How long should the rasters be in milliseconds.
    :return:
    """

    if isinstance(clu, str):
        clu = h5.get_node(clu)

    fs = h5.root.events._v_attrs['sample_rate_Hz']
    n_starts = len(start_sample_array)
    len_array = np.int(time_window_ms * 0.001 * fs)
    rasters = np.empty((n_starts, len_array))
    rasters[:] = np.nan
    clu = clu[:]
    time_window_samples = time_window_ms * 0.001 * fs
    for i, start in enumerate(start_sample_array):
        end = np.int(time_window_samples + start)
        a = clu[:] >= start
        b = clu[:] <= end
        c = a * b
        n = np.sum(c)
        ras = clu[c]
        rasters[i, :n] = ras - start
    if convert_to_ms:
        return rasters / (fs / 1000.)
    else:
        return rasters



# @numba.autojit('i8(f8[:,:], f8, f8')
@numba.jit
def make_spike_array(raster_array, binsize, start_time, end_time):
    bin_edges = np.arange(start_time, end_time, binsize)  # will automatically truncate if time window is not divisible by binsize.
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    spike_array = np.zeros((len(bin_edges)-1, len(raster_array)), dtype=np.int64)
    # print bin_edges
    for i in xrange(len(raster_array)):
        raster = raster_array[i, :]
        raster = raster[~np.isnan(raster)]
        for ii in xrange(len(bin_edges)-1):
            bin_start = bin_edges[ii]
            bin_end = bin_edges[ii+1]
            raster_l = (raster >= bin_start) * (raster < bin_end)
            spike_array[ii, i] = np.sum(raster_l)
    return spike_array, bin_centers