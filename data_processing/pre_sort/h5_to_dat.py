#!/home/chris/anaconda/bin/python2.7
from os import path

__author__ = 'chris'


import tables
import numpy as np
from utils.probe_utils import load_probe
import sys


def h5_to_bin(raw_kwd_file_path, probe_file=None):
    """

    :param raw_kwd_file_path:
    :param probe_file:
    :return:
    """
    nchannels = 64
    if probe_file:
        prb = load_probe(probe_file)
        keys = prb.keys()
        keys.sort()
        ch_map = np.array([])
        for k in keys:
            grp = prb[k]
            # print grp
            graph = grp['graph']
            chs = []
            for i in graph:
                for ii in i:
                    chs.append(ii)
            cha = np.unique(chs)  # this should be in some kind of order if the graph was done reasonably..
            e = ch_map.size + cha.size
            print('Mapping in bin for shank %i:' %k)
            print('start ch: %i, end ch: %i'% (ch_map.size, e))
            ch_map = np.hstack((ch_map, cha))
    elif nchannels:  # assume first x channels.
        a = range(0, nchannels)
        ch_map = np.array(a)

    else:  # probe_file is None and nchannels=0:
        raise ValueError('number of channels must be specified by probe or nchannels argument')

    # load kwd:
    f = tables.open_file(raw_kwd_file_path, 'r')
    fn = path.splitext(raw_kwd_file_path)
    recordings = sorted([n._v_name for n in f.list_nodes('/recordings')])

    # construct array:
    for rec in recordings:
        print rec
        print type(rec)
        n = f.get_node('/recordings/%s/data'%rec)
        nrows = n.shape[0]
        #allocate huge memory:
        if probe_file is None:  # just make a range from 0 to nchannels in the recording if no probe.
            nchannels = n.shape[1]
            a = range(0, nchannels)
            ch_map = np.array(a)
        reconstructed = np.zeros((nrows, ch_map.size), dtype=np.int16)

        for ii, ch_i in enumerate(ch_map):
            print 'Copying ch %i' % ii
            reconstructed[:, ii] = n[:, ch_i]  #read a row into the correct spot in reconstructed matrix.

        # save array as flattened binary file:
        fn = path.splitext(raw_kwd_file_path)[0]
        save_fn = fn + '.%s' %rec
        save_fn = save_fn + '.dat'
        print 'Saving as: %s...' % save_fn
        reconstructed.tofile(save_fn, sep='')
        print 'Complete!\n'


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2:
        hlp_str = ("\nUnpacks *.raw.kwd to binary file (.dat) for viewing in neuroscope. "
                   "\n\nUsage: python h5_to_dat.py kwd_filename prb_filename."
                   "\nProbe file is used to order channels by shank.")
        sys.exit(hlp_str)
    else:
        raw_fn = args[0]
        prb_fn = args[1]

        h5_to_bin(raw_fn, prb_fn)
    sys.exit()