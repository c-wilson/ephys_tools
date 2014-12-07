__author__ = 'chris'

import copy
import param_util
import os

# acquisition_system = ''

# acquisition_system=''

def remap_sites_to_channels(site_graph, chan_translation):
    """
    Translate the site number to the recording channel number. Uses ordered sites list and corresponding ordered sgl_chan list.

    :param site_graph: list of lists of site-pairs (acquisition system channels are base 0 indexes!)
    :param chan_translation: dict of lists: {'sites': [], 'sgl_chan': []}.

    """

    # make translation for acquisition system (should be defined as global).
    if acquisition_system == '':
        ct = copy.deepcopy(chan_translation['sites'])
        # NOTE: assumes that site maps start with site 1, not with site 0!!!. Also that channels start at 0, not 1,
        # which is a good assumption if you aren't using matlab.
        trans = [x-1 for x in ct]
    else:
        trans = chan_translation[acquisition_system]

    chan_graph = copy.deepcopy(site_graph)  #don't want to change original site graph object.
    for i, pair in enumerate(site_graph):
        for ii, site in enumerate(pair):
            idx = chan_translation['sites'].index(site)
            chan_graph[i][ii] = trans[idx]
    return chan_graph


def calc_channel_list(site_graph, chan_translation, bad_channels=[]):
    """
    calculates the channel which a site was recorded on, and makes a 1d list of the unique channels in the map.
    Does not add channels enumerated in the bad channels parameter to the final channels list.

    :param site_graph: list of lists of site-spatial pairs.
    :param bad_sites: list of channels
    :return: list of channels
    """
    # TODO: remove bad channels from graph (klusta crashes when more channels in graph than in chan list).
    chan_graph = remap_sites_to_channels(site_graph, chan_translation)
    chan_list = []
    for pair in chan_graph:
        for chan in pair:
            if chan not in chan_list and chan not in bad_channels:
                chan_list.append(chan)
    chan_list.sort()  # sort in place.
    return chan_list


def make_geometry_by_channels(geo_by_site, chan_translation):
    """

    :param geo_by_site:
    :param chan_translation:
    :return:
    """

    # make translation for acquisition system (should be defined as global).
    if acquisition_system == '':
        ct = copy.deepcopy(chan_translation['sites'])
        # NOTE: assumes that site maps start with site 1, not with site 0!!!.
        trans = [x-1 for x in ct]
    else:
        trans = chan_translation[acquisition_system]

    sites = geo_by_site.keys()
    geo_by_chan = {}
    for site in sites:
        idx = chan_translation['sites'].index(site)
        ch = trans[idx]
        geo_by_chan[ch] = geo_by_site[site]
    return geo_by_chan


def load_probe(filename, acq_system=''):
    """
    loads channel groups in probe file as a variable. uses a 'global' namespace unlike kk, so you can use functions.
    :param filename: fullfilename to prb file.
    :return: channel_groups dict.
    """
    global acquisition_system
    acquisition_system = acq_system
    prb = {}
    execfile(filename, globals(), prb)

    return prb['channel_groups']


def write_probe(prb, filename):
    """
    writes a probe file adhering to conventions of kk. (removes functions).

    :param prb: probe 'channel_groups' object
    :param filename: filename to write (str)
    :return: none
    """
    orig_fn = filename + '.orig'
    if not os.path.isfile(orig_fn) and os.path.isfile(filename):
        os.rename(filename, orig_fn)
    prbdict = {'channel_groups': prb}
    prbstr = param_util.pydict_to_python(prbdict)
    prbstr = '""" PROBE FILE GENERATED BY PROCESS MODULE """\n\n'+ prbstr
    f = open(filename, 'w')
    f.write(prbstr)
    f.close()
    return


def translate_prb(filename, outfile, acq_system=''):
    """

    :param filename: file to translate
    :param outfile: file out path
    :param acq_system: System on which data were acquired (matches key in sites to channels dictionary).
    :return:
    """
    prb = load_probe(filename, acq_system)
    write_probe(prb, outfile)
