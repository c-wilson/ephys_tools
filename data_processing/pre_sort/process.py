#!/home/chris/anaconda/bin/python2.7
# ^ this makes this file executable.
__author__ = 'chris'

import os
from copy import deepcopy
import time
import argparse
import numpy as np
import tables
import \
    spikedetekt2.core.script as klusta_entry  # this is important and may change with subsequent versions of the klustering suite!
from utils import param_util
import utils.probe_util as prb_utils
import sys




class PreProcessSess(object):
    """
    Classdocs
    """


    def __init__(self, sess_prm_ffn, force_reprocess=False):
        """

        :param sess_prm_ffn: fullfilename string for session .prm file
        :param force_reprocess: forces recreation of KWD files if they already exist for a given rec (Not yet implemented)
        :return:
        """
        self.path = os.path.split(sess_prm_ffn)[0]
        self.prm_ffn = sess_prm_ffn  # this is a fullfilename.
        self.parameters = param_util.get_params(sess_prm_ffn)
        self.prb = prb_utils.load_probe(os.path.join(self.path, self.parameters['prb_file']),
                                        acq_system=self.parameters['acquisition_system'])
        print prb_utils.acquisition_system
        prb_utils.write_probe(self.prb, os.path.join(self.path, self.parameters['prb_file']))
        self.rec_prm_fns = {}
        self.generate_rec_prms()
        self.recs = {}
        self.process_recs()
        # TODO: make it possible to reprocess from already generated raw.kwd files without rewriting them!
        return

    def generate_rec_prms(self):
        for rec, runs in self.parameters['raw_data_files'].iteritems():
            tmp_rec_prms = deepcopy(self.parameters)
            tmp_rec_prms['raw_data_files'] = runs
            tmp_rec_prms['experiment_name'] = (self.parameters['experiment_name'] + '_rec_' + str(rec))
            tmp_fn = tmp_rec_prms['experiment_name'] + '.prm'
            tmp_ffn = os.path.join(self.path, tmp_fn)
            self.rec_prm_fns[rec] = tmp_ffn
            # TODO: add handling for recs that have individual parameters. Individual recs can be re-processed later if desire
            f = open(tmp_ffn, 'w')  # make file with overwrite.
            f.write(param_util.pydict_to_python(tmp_rec_prms))
            f.close()
        return

    def process_recs(self):
        for rec, fn in self.rec_prm_fns.iteritems():
            # TODO: think about parallelizing this?? (complicated: needs to factor ncores available and memory available-(nshanks*nrecs))
            stime = time.time()
            print 'Preprocessing rec ' + str(rec)
            self.recs[rec] = PreProcessRec(fn)
            etime = time.time() - stime
            print 'Rec ' + str(rec) + ' processing completed in ' + str(etime) + ' sec.'


            # try:
            #     stime = time.time()
            #     print 'Preprocessing rec ' + str(rec)
            #     self.recs[rec] = PreProcessRec(fn)
            #     etime = time.time() - stime
            #     print 'Rec '+ str(rec) + ' processing completed in ' + str(etime) + ' sec.'
            # except Exception as msg:
            #     etime = time.time() - stime
            #     print '*ERROR* Rec '+ str(rec) + ' processing failed in ' + str(etime) + ' sec.'
            #     print msg
        print 'Session processing completed.'
        return

    def run_klusta(self, **kwargs):
        """
        runs run_klusta() for each rec instance.
        :param kwargs: to pass to su
        :return: none
        """
        for key, rec in self.recs.iteritems():
            try:
                print 'Running klusta for rec ' + str(key)
                rec.run_klusta(**kwargs)
            except Exception as msg:
                print 'ERROR klusta run failed for rec ' + str(key)
                print msg
        return


class PreProcessRec(object):
    """
    Classdocs
    """
    #TODO: copy voyeur edata here.
    def __init__(self, prm_ffn, **kwargs):
        """

        :param prm_ffn: fullfilename to .prm file for the rec.
        :param kwargs: not implemented
        :return:
        """
        self.path = os.path.split(prm_ffn)[0]
        self.prm_ffn = prm_ffn  # this is a fullfilename!
        self.parameters = param_util.get_params(filename=prm_ffn)
        self.prb = prb_utils.load_probe(os.path.join(self.path, self.parameters['prb_file']),
                                        acq_system=self.parameters['acquisition_system'])
        prb_utils.write_probe(self.prb, os.path.join(self.path, self.parameters['prb_file']))
        # TODO: reuse raw.kwd files without rewriting.
        self.data_ffn = self.path + self.parameters['experiment_name'] + '.raw.kwd'
        self.data_file = tables.open_file(self.data_ffn, 'w')  # CURRENTLY OVERWRITES!!!
        self.data_file.create_group('/', 'recordings')
        self.run_ephys_fns = []
        self.run_beh_fns = []
        for run in self.parameters['raw_data_files']:
            self.run_ephys_fns.append(run[0])
            self.run_beh_fns.append(run[1])
        self.runs = []
        self.append_runs()
        self.data_file.close()
        self.update_prm_file()
        return

    def update_prm_file(self):
        """
        updates prm file to use the new filename.
        :return:
        """
        self.parameters['raw_data_files'] = self.data_ffn
        prms_str = param_util.pydict_to_python(self.parameters)
        f = open(self.prm_ffn, 'w')
        f.write(prms_str)
        f.close()
        return

    def append_runs(self):
        """
        this is the meat of the matter.
        :return:
        """
        if isinstance(self.run_ephys_fns, str):
            self.run_ephys_fns = [self.run_ephys_fns]  # make into list
        for i, (run_ephys_fn, run_beh_fn) in enumerate(zip(self.run_ephys_fns, self.run_beh_fns)):
            # run_data_ffn = os.path.join(self.rec_path, run_data_fn)
            # print '\tPreprocessing run ' + str(i)
            # run_grp = self.data_file.create_group('/recordings',str(i))
            # self.runs.append(PreProcessRun(run_data_ffn, run_grp, self.data_file, self.prms, self.prb))
            run_ephys_fn = os.path.join(self.path, run_ephys_fn)
            run_beh_fn = os.path.join(self.path, run_beh_fn)
            print '\tPreprocessing run ' + str(i)
            run_grp = self.data_file.create_group('/recordings', str(i))
            self.runs.append(PreProcessRun(run_ephys_fn, run_beh_fn, run_grp, self.data_file, self.parameters, self.prb))
            # try:
            #     # run_grp = self.data_file.create_group('/recordings', str(i))
            #
            #
            #     # self.runs.append(PreProcessRun(run_data_ffn, run_grp, self.data_file, self.parameters, self.prb))
            # except Exception as msg:
            #     print 'ERROR: problem with run: ' + run_data_fn
            #     print msg


    def run_klusta(self, **kwargs):
        """

        :param kwargs:
        :return:
        """
        #TODO: add this back with subprocess support.
        # klusta_entry.run_all(self.prm_ffn, **kwargs)
        pass


class PreProcessRun(object):
    def __init__(self, run_bin_fn, run_beh_fn, run_grp, rec_h5_obj, rec_prms, prb):
        """

        :param run_bin_fn:
        :param run_beh_fn:
        :param run_grp:
        :param rec_h5_obj:
        :param rec_prms:
        :param prb:
        :return:
        """
        # TODO: move the pl_trig_chan calculation to the rec class.
        assert isinstance(rec_h5_obj, tables.File)
        self.beh_fn = run_beh_fn
        self.prms = rec_prms
        self.run_group = run_grp
        self.bin_fn = run_bin_fn
        meta_fn = self.bin_fn[:-3] + 'meta'
        self.sgl_meta = read_meta_file(meta_fn)
        self.channels, self.chan_idxes = calc_channels(self.sgl_meta, rec_prms)  # shape = (channels,)
        self.nchannels = len(self.channels)
        self.ephys_channels = calc_ephys_channels(prb)
        nephys_channels = len(self.ephys_channels)
        self.rec_h5_obj = rec_h5_obj
        if rec_prms['nchannels'] != nephys_channels:  #nchannels in rec parameters file is only neural!
            print ('WARNING: Number of neural channels specified in .prm file does not match the number of ' \
                   'channels specified in the .prb file.')
        # self.add_voyeur_behavior_data()
        self.edata = self.add_raw_ephys_data()

        if 'pl_trigger' in self.edata.keys():
            self.data_plfilt = self.rm_pl_noise(rec_h5_obj, rec_prms['pl_trig_chan'])
        else:
            print 'NO powerline trigger signal found, using raw data without filtering.'
            self.edata['neural'].rename('data')  # make the raw data into the data stream.
            self.data_plfilt = self.edata

    def add_raw_ephys_data(self, max_load=1e9):
        """

        :param rec_h5_obj: h5 file object from the record.
        :param max_load: number of integers to load at a given time for RAM limitations (default uses 8 GB).
        :return:
        """

        filesize = os.path.getsize(self.bin_fn)
        expct_rows = filesize / self.nchannels / 2
        data = {}  # dictionary to hold all of the data array objects (neural and metadata streams).
        for k, v in self.chan_idxes.iteritems():
            data[k] = self.rec_h5_obj.create_earray(self.run_group,
                                                    name=k,
                                                    atom=tables.Int16Atom(),
                                                    shape=(0, len(v)),
                                                    title='raw %s edata' % k,
                                                    expectedrows=expct_rows)
            data[k]._v_attrs['bin_filename'] = self.bin_fn

        f = open(self.bin_fn, 'rb')
        ld_q = int(max_load) / int(self.nchannels)  # automatically floors this value. a ceil wouldn't be bad
        ld_iter = ld_q * self.nchannels  # calculate number of values to read in each iteration
        ld_count = 0
        print '\t\tAdding raw run recording data to kwd...'
        while ld_count < filesize:
            arr = np.fromfile(f, np.int16, ld_iter)
            ld_count += ld_iter
            larr = arr.size / self.nchannels
            arr.shape = (larr, self.nchannels)
            for k, v in data.iteritems():
                idx = self.chan_idxes[k]
                v.append(arr[:, idx])
                v.flush()
            pc = float(ld_count)/float(filesize) * 100.
            if pc > 100.:
                pc = 100.
            print '\t\t\t... %0.1d %% complete' % pc
        f.close()
        self.run_group._v_attrs['bin_filename'] = str(self.bin_fn)
        self.rec_h5_obj.flush()
        return data

    def add_voyeur_behavior_data(self):
        """
        Copies behavior data from Voyeur HDF5 into run h5 group.

        :return:
        """
        print '\t\tAdding Voyeur data to kwd...'
        beh_h5 = tables.open_file(self.beh_fn, mode='r')
        # will copy the voyeur data along with metadata attributes.
        beh_group = self.rec_h5_obj.copy_node(beh_h5.root, self.run_group, newname='Voyeur_data', recursive=True)
        self.run_group._v_attrs['Voyeur_filename'] = str(self.beh_fn)  # save the filename for later, just in case...
        self.rec_h5_obj.flush()

        beh_h5.close()


    def rm_pl_noise(self, rec_h5_obj, pl_trig_chan):
        """

        :param rec_h5_obj:
        :param pl_trig_chan: SpikeGL channel  number where the PL trigger is recorded. THIS IS SpikeGL channel,
                             not index in the edata matrix!
        :return: pointer to H5 array.
        """
        data_raw = self.edata['neural']
        filtered_data_array = rec_h5_obj.create_earray(self.run_group, 'data',
                                                  atom=tables.Int16Atom(),
                                                  shape=(data_raw.shape[0], 0),
                                                  title='pl trigger (60 Hz) filtered neural data',
                                                  expectedrows=data_raw.shape[1])
        print '\t\tPL filtering neural channels. Loading pl_trigger.'
        pl_trig_sig = self.edata['pl_trigger'].read()[:, 0]  # array is multidimensional, so we just want the first column
        threshold = np.mean(pl_trig_sig)
        pl_trig_log = pl_trig_sig > threshold
        pl_edge_detect = np.convolve([1, -1], pl_trig_log, mode='same')
        pl_edge_idx = np.where(pl_edge_detect[1:] == 1)[0]  # find upward edge - np.where returns tuple, want 1st dim
        pl_len = np.diff(pl_edge_idx).max()
        n_pl = pl_edge_idx.size
        sig_pl = np.zeros((n_pl - 1, pl_len), dtype=np.int16)
        n_ch = data_raw.shape[1]
        # Filter and save every channel in a loop. Since only one channel is loaded and processed at a time, this will
        # use only as much memory as loading a single channel of data at a time at the expense of speed.
        for ch_i in xrange(n_ch):
            chan_sig = data_raw[:, ch_i]
            print '\t\t\tfiltering channel %i of %i...' % (ch_i, n_ch)
            sig_len = chan_sig.size
            chan_sig = chan_sig - chan_sig.mean()
            for i, edge in enumerate(pl_edge_idx):
                end = edge + pl_len
                if end < sig_len:
                    sig_pl[i, :] = chan_sig[edge:end]
            pl_template = sig_pl.mean(axis=0)
            for i in xrange(n_pl - 1):
                st = pl_edge_idx[i]
                end = pl_edge_idx[i + 1]
                l = end - st
                chan_sig[st:end] -= pl_template[:l]
            filtered_data_array.append(chan_sig[:, np.newaxis])  # save processed once complete with every channel.
            filtered_data_array.flush()
        return filtered_data_array



def calc_ephys_channels(prb):
    """

    :param prb: probe object, passed from session or rec.
    :return: list of ephys channels
    """
    ephys_channels = np.array([], dtype=int)
    for shank in prb.values():
        ephys_channels = np.append(ephys_channels, shank['channels'])
    ephys_channels.sort()
    return ephys_channels


def calc_channels(sgl_meta, prms):
    """

    :param sgl_meta: sgl_meta dictionary from read_meta_file()
    :return:
    """

    chanstr = sgl_meta['saveChannelSubset']
    chanstr_sp = chanstr.split(',')
    chan_arr = np.array([], dtype=np.int)
    for st in chanstr_sp:
        if not st.find(':') == -1:
            start = int(st.split(':')[0])
            stop = int(st.split(':')[1])
            chans = np.arange(start, stop + 1)
            chan_arr = np.append(chans, chan_arr)
        else:
            chan_arr = np.append(chan_arr, int(st))
    channel_to_idx = {}
    for k, v in prms['chan_config'].iteritems():
        if isinstance(v, int):  #allows us to iterate.
            v = [v]
        trans = []
        try:
            for ch in v:
                    trans.append(np.where(chan_arr==ch)[0][0])
            channel_to_idx[k] = trans
        except IndexError:  # this will happen when the channel is not found in the channel array (from meta file)
            print('WARNING: channel %i was not recorded, according to spikeGL meta file! '
                  'No "%s" signal will be copied' % (ch, k))

    return chan_arr, channel_to_idx


def read_meta_file(fn):
    """
    For reading files adhering to spikeGL sgl_meta file output format.

    :param fn: filename.sgl_meta
    :return: dictionary containing sgl_meta edata.
    """

    f = open(fn, 'r')
    st = f.read()
    st_sp = st.splitlines()
    out = {}
    for ln in st_sp:
        lnsp = ln.split(' = ')
        out[lnsp[0]] = lnsp[1]
    return out


def main():
    """
    Entry point for chris's preprocessing.
    :return:
    """
    # TODO: send output to log file (or add the option to do so).
    dtg = time.strftime('%Y%m%d_%H%M%S')
    sys.stdout = open('process_log_%s.log' % dtg, 'w', buffering=10)

    if not klusta_entry.check_path():  # check that klustakwik is in path.
        return

    parser = argparse.ArgumentParser(description='Run preprocessing for kk3.')
    parser.add_argument('prm_filename',
                        help='.prm filename')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='only run clustering on first few seconds of edata')
    parser.add_argument('--overwrite', action='store_true', default=False,
                        help='overwrite existing KWIK files if they exist')
    parser.add_argument('--detect-only', action='store_true', default=False,
                        help='run only spikedetekt')
    parser.add_argument('--cluster-only', action='store_true', default=False,
                        help='run only klustakwik (after spikedetekt has run)')
    args = parser.parse_args()
    runsd, runkk = True, True
    if args.detect_only:
        runkk = False
    if args.cluster_only:
        runsd = False

    prms = param_util.get_params(args.prm_filename)
    if type(prms['raw_data_files']) is dict:
        pp = PreProcessSess(args.prm_filename)
    else:
        pp = PreProcessRec(args.prm_filename)

    print '\n\nPREPROCESS COMPLETE, running klustas'


    return None


if __name__ == "__main__":
    main()