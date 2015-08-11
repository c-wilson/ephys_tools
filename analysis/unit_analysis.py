from __future__ import division
import tables as tb
import numpy as np
import matplotlib as plt
from .data_classes import Recording
from utils.h5_decorator import h5decorator
import numba as nb

# import numba as b

__author__ = 'chris'

@h5decorator
def plot_spikes(h5, clu):
    pass


@h5decorator
def calc_autocorrelelogram(h5, clu, half_width_ms=20):
    """
    Calculates an autocorrelelogram for a given cluster. Returns an array of size half_width_ms*2 which is symmetric
    around the half_width_ms index:

        ac[:half_width_ms] == np.flipud(ac[half_width_ms])

    Can be nicely plotted with: plt.step(np.arange(-half_width_ms , half_width_ms), ac/ac.sum(), where ='post')

    :param h5:
    :param clu:
    :param half_width_ms:
    :return: np.int array of autocorrelelogram
    """

    if isinstance(clu, str):
        clu = h5.get_node(clu)
    ac = np.zeros(half_width_ms*2, dtype=np.int)
    clu_ms = sample_to_ms(h5, clu.read())  #returns float
    _acorr2(ac, clu_ms, int(half_width_ms))  # edits the ac array in place and returns none.

    # make the first half of the ac (symmetry is defined for ac):
    for i in xrange(half_width_ms):
        ac[i] = ac[-(1+i)]

    return ac


@nb.jit(nb.void(nb.int64[:], nb.float64[:], nb.int64), nopython=True)
def _acorr2(corr_array, clu_ms, half_width):
    """
    Private function to allow numba nopython. Calculates half of the autocorrelelogram. Modifies an array of zeros in
    place. This is massively fast in numba even if it is more optimized for pure python (~81 ms vs 373 ns).

    :param corr_array: THIS IS MODIFIED INPLACE
    :param clu_ms:
    :param half_width:
    :return: Nothing!
    """
    nspikes = len(clu_ms)
    for i in xrange(nspikes):
        spikei = clu_ms[i]
        for j in xrange(i+1, nspikes):
            spikej = clu_ms[j]
            diff = spikej - spikei
            if diff < half_width:
                i_diff = np.int(np.floor(diff))
                i_diff += half_width
                corr_array[i_diff] += 1
            else:
                break  # breaks to outer loop.
    return




def sample_to_ms(h5, samples, sample_rate=20833.):
    """
    Convert sample time value(s) to millisecond timebase.

    :param samples: numeric (array or scalar)
    """

    sample_rate = h5.root.events._v_attrs['sample_rate_Hz']
    ms = samples/sample_rate * 1000

    return ms