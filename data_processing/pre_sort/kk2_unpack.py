__author__ = 'chris'

import sys
import os.path
import tables
import numpy as np


def unpack_kk2(kwik_path, kwx_path, output_dir):
    """
    Unpacks klustakwik/spikedetect2 h5 files (kwik and kwd) files into legacy KK/NDManager filetypes for use with
    legacy manual klustering programs.

    :param argv: tuple containing ('kwik_file_path.kwik', 'kwx_file_path.kwx', 'output_files_dir')
    :return:
    """

    filenamebase = os.path.splitext(os.path.split(kwik_path)[1])[0]
    output_base_path = os.path.join(output_dir, filenamebase)  # output filename without extension.
    kwik = tables.open_file(kwik_path, mode='r')
    kwx = tables.open_file(kwx_path, mode='r')
    res_dict = {}
    filtered_h5_path = os.path.join(output_dir, filenamebase) + '.high.kwd'
    
    for group in kwik.root.channel_groups:
        grp_num = int(group._v_name) + 1  # old system expects that shanks start at 1, not 0.
        clu_out_path = output_base_path + '.clu.' + str(grp_num)
        res_out_path = output_base_path + '.res.' + str(grp_num)

        # Handle clu.
        clu_array = group.spikes.clusters.original.read()
        n_clusters = clu_array.max() + 1  # first cluster is '0'.
        clu_file = open(clu_out_path, 'w')
        clu_file.write('%i\n' % n_clusters)
        np.savetxt(clu_file, np.int16(clu_array),fmt="%i")
        clu_file.close()

        # Handle res.
        res_array = group.spikes.time_samples.read()
        np.savetxt(res_out_path, res_array, fmt='%i')
        # we need to save spike times for later to add to the fet file array. CANNOT depend on order, so put in a dict.
        res_dict[grp_num] = res_array
    kwik.close()

    for group in kwx.root.channel_groups:
        grp_num = int(group._v_name) + 1  # old system expects that shanks start at 1, not 0.
        spk_out_path = output_base_path + '.spk.' + str(grp_num)
        fet_out_path = output_base_path + '.fet.' + str(grp_num)

        #  Handle spk.
        spk_array = group.waveforms_filtered.read()  # (spikes, samples, channels)
        # spk_array = spk_array.transpose(0,2,1)  # need (spikes, channels, samples) UPDATE: NO!!!
        spk_array.tofile(spk_out_path)

        # Handle fet.
        fet_node = group.features_masks
        if fet_node.ndim == 3 and fet_node.shape[2] == 2:
            fet_array = group.features_masks[:,:,0]
        if fet_array.dtype == np.float32:  # need to convert to int16 as per spikedetect algorithm.
            factor = 2.**12  # nonsensical number from spikedetect.
            factor = factor/np.abs(fet_array).max()
            if not (factor == 1):  # don't copy the array if we don't need to.
                fet_array = (fet_array * factor).astype(np.int16)
        res_tmp = res_dict[grp_num]  # spike times to add as the last column in the fet array.
        times_tmp = np.expand_dims(res_tmp, axis=1)  # times dtype is uint64.
        if np.max(times_tmp) >= 2**52:
            print 'WARNING PRECISION LOSS IN TIME FEATURE.'
        fet_array = np.concatenate((fet_array, times_tmp), axis=1)
        # WARNING: this results in an float64 array, may result in precision loss from the uint64 time expression if
        # recording is more than ~1250999896 hrs in duration (float64 has 52 bit significand).
        fet_file = open(fet_out_path, 'w')
        fet_file.write('%i\n' % fet_array.shape[1])  # first line is number of PCs
        np.savetxt(fet_file, fet_array, fmt="%i")  # add rows of array as lines of text.
        fet_file.close()
    kwx.close()

    # GENERATE FILTERED RAW .DAT FILE
    # if it doesnt exist skip this part
    if os.path.isfile(filtered_h5_path):
        filtered_h5 = tables.open_file(filtered_h5_path)
        nd = filtered_h5.get_node('/recordings/0/data')
        filt_arr = nd.read()  # results in an int16 samples x channels array.
        # save
        filter_out_path = output_base_path + '.fil'
        # print 'saving bin:'
        # print '%s' % filter_out_path
        filt_arr.tofile(filter_out_path, sep='') # sep='' makes binary file of dtype corresponding to array (here int16).
        # Data is flattened in C order by the 'tofile' function, which is what we want.
        # To read: np.fromfile(filename, dtype=np.int16)
        filtered_h5.close()
    else:
        print 'Hi-pass filtered file not found, skipping .fil file creation'

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 3:
        hlp_str = ("\nUnpacker for kwik/kwd files to legacy NDManager files (.fet, .clu, .spk, .res).\n" +
                   "Usage: python kk2_unpack.py 'kwik_file_path.kwik' 'kwx_file_path.kwx' 'output_files_dir'\n")
        sys.exit(hlp_str)
    else:
        kwk_p = args[0]
        kwx_p = args[1]
        out_dir = args[2]
        unpack_kk2(kwk_p, kwx_p, out_dir)
    sys.exit()
