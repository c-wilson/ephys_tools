__author__ = 'chris'

import copy


def remap_sites_to_channels(site_graph, chan_translation):
    """
    Translate the site number to the recording channel number. Uses ordered sites list and corresponding ordered sgl_chan list.

    :param site_graph: list of lists of site-pairs
    :param chan_translation: dict of lists: {'sites': [], 'sgl_chan': []}.

    """

    chan_graph = copy.deepcopy(site_graph)  #don't want to change original site graph object.
    for i, pair in enumerate(site_graph):
        for ii, site in enumerate(pair):
            idx = chan_translation['sites'].index(site)
            chan_graph[i][ii] = chan_translation['sgl_chan'][idx]
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
