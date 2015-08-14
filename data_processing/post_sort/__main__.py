__author__ = 'chris'
"""
Module combines post-sorted data from kwik and kwd sources to make a cohesive ephys unit.

"""

# TODO: need to propagate root metadata from the behavior file to the h5 file!!!

import tables as tb
import os
import argparse

from ephys_tools.data_processing.pre_sort.utils.param_util import get_params
import logging
from eventer import eventer
from clusterer import clusterer
from streamer import streamer
import warnings

warnings.simplefilter('ignore', tb.NaturalNameWarning)
logging.basicConfig(level=logging.INFO)

def main(input_filename, overwrite=False, append=False):
    """

    :param input_filename: KWIK or .prm file.
    :param overwrite: overwrite destination file if it exists.
    :param append: append changed clusters to existing file.
    :return:
    """
    ex = os.path.splitext(input_filename)[1]
    input_dir = os.path.split(input_filename)[0]

    if overwrite and append:
        raise ValueError(u'Post sorting can not be run with both --overwrite and --append flags')
    if ex.lower() == u'.prm':
        prms = get_params(input_filename)
        if isinstance(prms[u'raw_data_files'], dict):
            input_filebases = [prms[u'experiment_name'] + u'_rec_' + x for x in prms[u'raw_data_files']]
            input_kwiks = [os.path.join(input_dir, x) + u'.kwik' for x in input_filebases]
    elif ex.lower() == u'.kwik':
        #only a single kwik to deal with.
        # input_kwiks = [input_filename]
        t = os.path.split(input_filename)[1]
        input_filebases = [os.path.splitext(t)[0]]
    else:
        raise ValueError(u'Input filename was an unknown type: must be .kwik or .prm file.')


    # check that input files exist in the directory that we expect.
    filenames = []
    for fb in input_filebases:
        fbs = {}
        logging.debug(u'Checking for files with filebase {0:s}'.format(fb))
        for ext in (u'.raw.kwd', u'.kwik'):
            file = os.path.join(input_dir, fb) + ext
            fbs[ext] = file
            if not os.path.exists(file):
                raise FileException(u'File does not exist: {0:s}'.format(file))
            else:
                logging.debug(u'File {0:s} exists!'.format(file))
        fbs[u'destination'] = os.path.join(input_dir, fb) + u'.h5'
        logging.info(u'All {0:s} input files present.'.format(fb))
        filenames.append(fbs)


    # check that files don't exist or that we're either appending or overwriting them. Do this first, so that we can
    # throw an exemption before we waste time processing any of the files.   for file in destination_filenames:
    for files in filenames:
        file = files[u'destination']
        if os.path.exists(file):
            logging.debug(u'Destination file exists: %s' % file)
            if not overwrite and not append:
                logging.error(
                    u'Destination file exists: {0:s} no overwrite or append parameters were passed aborting.'.format(
                        file))
                if len(filenames) > 1:
                    raise FileException(u'Destination file %s already exists. Please use --overwrite or --append to '
                                        u'reprocess or update file.')

    # if things are ok, put this thing together.
    for files in filenames:
        file = files[u'destination']
        if os.path.exists(file) and append:
            with tb.open_file(file, 'a') as dest_file:
                d_mod = dest_file.get_node_attr(u'/clusters', u'kwik_mod_time')
                if d_mod >= os.path.getmtime(files['.kwik']):
                    logging.info(u'Existing file contains most recent kwik data, skipping..'.format(file))
                else:
                    logging.info(u'Updating kwik data for {0:s}'.format(file))
                    clusterer(files['.kwik'], dest_file)
        elif (os.path.exists(file) and overwrite) or not os.path.exists(file):
            logging.info(u'Creating destination file: {0}'.format(file))
            with tb.open_file(file, 'w') as dest_file:
                clusterer(files['.kwik'], dest_file)
                eventer(files['.raw.kwd'], dest_file)
                streamer(files['.raw.kwd'], dest_file)
        else:
            raise Exception(u'File already exists and no overwrite parameter was provided:'
                            u'\n\tFile: {0}'
                            u'\n\tfile exists = {1}'
                            u'\n\tappend = {2}'
                            u'\n\toverwrite = {3}.'.format(file, os.path.exists(file), append, overwrite))


def append(destination_file, files):
    """

    :param destination_file:
    :param files:
    :return:
    """

def make(filebase, refresh=True):
    pass



def update(filebase):
    pass


class FileException(Exception):
    pass


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description='Run preprocessing for kk3.')

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

    main(args.input_filename, overwrite=args.overwrite, append=args.append)