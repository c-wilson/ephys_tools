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
from subprocess import Popen
import logging




class PreProcessSess(object):
    """
    Classdocs
    """


    def __init__(self, sess_prm_ffn, force_reprocess, klusta_args):
        """

        :param sess_prm_ffn: fullfilename string for session .prm file
        :param force_reprocess: forces recreation of KWD files if they already exist for a given rec (Not yet implemented)
        :return:
        """
        self.path = os.path.split(sess_prm_ffn)[0]
        # self.processed_path = os.path.join(self.path, '/processed')
        self.force_reprocess = force_reprocess
        self.prm_ffn = sess_prm_ffn  # this is a fullfilename.
        self.parameters = param_util.get_params(sess_prm_ffn)
        self.prb = prb_utils.load_probe(os.path.join(self.path, self.parameters['prb_file']),
                                        acq_system=self.parameters['acquisition_system'])
        # print prb_utils.acquisition_system
        prb_utils.write_probe(self.prb, os.path.join(self.path, self.parameters['prb_file']))
        self.rec_prm_fns = {}
        self.generate_rec_prms()
        self.recs = {}
        self.process_recs(klusta_args)
        # TODO: make system move raw data files to new directory or make process files in new directory.
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

    def process_recs(self, klusta_args):
        for rec, fn in self.rec_prm_fns.iteritems():
            stime = time.time()
            logging.info('Preprocessing rec ' + str(rec))
            try:
                self.recs[rec] = PreProcessRec(fn, self.force_reprocess, klusta_args)
                etime = time.time() - stime
                logging.info('Rec ' + str(rec) + ' preprocessing completed in ' + str(etime) + ' sec.')
            except Exception as msg:
                logging.exception('error processing rec %s:'%rec)
                print 'error processing rec %s'%rec
        logging.info('Session processing completed.')
        return


class PreProcessRec(object):
    """
    Classdocs
    """
    def __init__(self, prm_ffn, force_reprocess, klusta_args, **kwargs):
        """

        :param prm_ffn: fullfilename to .prm file for the rec.
        :param kwargs: not implemented
        :return:
        """
        self.klusta_instance = None
        self.path = os.path.split(prm_ffn)[0]
        self.force_reprocess = force_reprocess
        self.prm_ffn = prm_ffn  # this is a fullfilename!
        self.parameters = param_util.get_params(filename=prm_ffn)
        self.prb = prb_utils.load_probe(os.path.join(self.path, self.parameters['prb_file']),
                                        acq_system=self.parameters['acquisition_system'])
        prb_utils.write_probe(self.prb, os.path.join(self.path, self.parameters['prb_file']))
        self.data_ffn = self.path + self.parameters['experiment_name'] + '.raw.kwd'
        if self.force_reprocess or not os.path.exists(self.data_ffn):
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
        #TODO: add case for not checking if PL filtered data is contained in H5 and re-making it from raw if not.
        self.update_prm_file()
        self.run_klusta(klusta_args)
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
            run_ephys_fn = os.path.join(self.path, run_ephys_fn)
            run_beh_fn = os.path.join(self.path, run_beh_fn)
            logging.info( '\tPreprocessing run ' + str(i))
            run_grp = self.data_file.create_group('/recordings', str(i))
            self.runs.append(PreProcessRun(run_ephys_fn, run_beh_fn, run_grp, self.data_file, self.parameters, self.prb))


    def run_klusta(self, klusta_args, **kwargs):
        """

        :param kwargs:
        :return:
        """

        if not (klusta_args['runkk'] + klusta_args['runsd']):
            loggind.info('Skipping klusta for %s.' %self.data_ffn)
            return

        error_out_fn = '%s_klustalog.log'%self.data_ffn
        cmd = 'klusta %s' %  self.prm_ffn
        if klusta_args['runkk'] and not klusta_args['runsd']:
            cmd += ' --cluster-only'
        elif klusta_args['runsd'] and not klusta_args['runkk']:
            cmd += ' --detect-only'

        if klusta_args['overwrite']:
            cmd += ' --overwrite'

        logging.info('opening subprocess: %s' %cmd)
        with open(os.devnull, 'w') as devnull:
            self.klusta_instance = Popen(cmd,
                                         shell=True,
                                         stdout=devnull,
                                         stderr=open(error_out_fn, 'w', buffering=0))
        self.klusta_pid = self.klusta_instance.pid
        logging.info('\n ** Running klusta on %s. Running on PID: %i' %(self.data_ffn, self.klusta_pid))

        return

    def wait_klusta(self):
        if self.klusta_instance is not None:
            exit_code = self.klusta_instance.wait()
            if exit_code != 0:
                logging.error('Klusta instance for %s completed with exit code %i' % (self.data_ffn, exit_code))
            else:
                logging.info('Klusta instance for %s completed with exit code %i' % (self.data_ffn, exit_code))
                #TODO: delete PL filtered data node from run groups
        else:
            logging.warning('No klusta instance initiated for %s' % self.data_ffn)



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
            logging.warning('Number of neural channels specified in .prm file does not match the number of ' \
                   'channels specified in the .prb file.')
        self.add_voyeur_behavior_data()
        self.edata = self.add_raw_ephys_data()
        if 'pl_trigger' in self.edata.keys():
            self.data_plfilt = self.rm_pl_noise(rec_h5_obj, rec_prms['pl_trig_chan'])
        else:
            logging.warning('NO powerline trigger signal found, using raw data without filtering.')
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
            data[k]._v_attrs['acquisition_system'] = self.prms['acquisition_system']
            data[k]._v_attrs['sampling_rate_Hz'] = self.prms['sample_rate']

        f = open(self.bin_fn, 'rb')
        ld_q = int(max_load) / int(self.nchannels)  # automatically floors this value. a ceil wouldn't be bad
        ld_iter = ld_q * self.nchannels  # calculate number of values to read in each iteration
        ld_count = 0
        logging.info('\t\tAdding raw run recording data to kwd...')
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
            logging.info('\t\t\t... %0.1d %% complete' % pc)
        f.close()
        self.run_group._v_attrs['bin_filename'] = str(self.bin_fn)
        self.rec_h5_obj.flush()
        return data

    def add_voyeur_behavior_data(self):
        """
        Copies behavior data from Voyeur HDF5 into run h5 group.

        :return:
        """
        logging.info('\t\tAdding Voyeur data from %s to kwd...'%self.beh_fn )
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
        logging.info('\t\tPL filtering neural channels. Loading pl_trigger.')
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
            logging.info('\t\t\tfiltering channel %i of %i...' % (ch_i, n_ch))
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
            logging.warning('WARNING: channel %i was not recorded, according to spikeGL meta file! '
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

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def print_path():
    print '\n'.join(os.environ["PATH"].split(os.pathsep))

def check_path():
    prog = 'klusta'
    if not (which(prog) or which(prog + '.exe')):
        print("Error: '{0:s}' is not in your system PATH".format(prog))
        return False
    return True


def main():
    """
    Entry point for chris's preprocessing.
    :return:
    """
    # TODO: send output to log file (or add the option to do so).

    # try:  # check that klusta is accessible before continuing.
    #
    #     t = Popen('klusta', stderr=open(os.devnull), shell=True)
    #     t.wait()
    # except OSError:
    #     raise EnvironmentError('klusta executable not found, make sure you are in the klusta environment before running!')

    if not check_path():  # check that klustakwik is in path.
        return

    parser = argparse.ArgumentParser(description='Run preprocessing for kk3.')
    parser.add_argument('prm_filename',
                        help='.prm filename')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='only run clustering on first few seconds of edata')
    parser.add_argument('--repreprocess', action='store_true', default=False,
                        help='overwrite existing raw.kwd file if they exist.')
    parser.add_argument('--preprocess_only',action='store_true', default=False,
                        help='only run preprocess routine, not klusta')
    parser.add_argument('--detect-only', action='store_true', default=False,
                        help='run only spikedetekt')
    parser.add_argument('--cluster-only', action='store_true', default=False,
                        help='run only klustakwik (after spikedetekt has run)')
    parser.add_argument('--overwrite', action='store_true', default=False,
                       help='overwrite the KWIK files is they already exist')
    args = parser.parse_args()
    runsd, runkk = True, True
    if args.detect_only:
        runkk = False
    if args.cluster_only:
        runsd = False
    if args.preprocess_only:
        runsd = False
        runkk = False

    klusta_args = {'runkk': runkk, 'runsd': runsd, 'overwrite': args.overwrite}

    prms = param_util.get_params(args.prm_filename)
    if type(prms['raw_data_files']) is dict:
        pp = PreProcessSess(args.prm_filename, args.repreprocess, klusta_args)

    res = [x.wait_klusta() for x in pp.recs.values()]



    return None


if __name__ == "__main__":
    dtg = time.strftime('%Y%m%d_%H%M%S')
    logging.basicConfig(filename='process_%s.log'%dtg, level=logging.INFO)

    main()

