__author__ = 'chris'

"""
System to get clusters into ephys file.

"""

import tables as tb
import time
import os.path
import logging
import numpy as np


def clusterer(kwik_fn, destination):
    """

    :param kwik_fn:
    :param destination: tables.File object that will be the destination for the cluster data from the .kwik file.
    :return:
    """
    logging.info('adding clusters')
    assert isinstance(destination, tb.File)
    dtg = time.strftime('D%Y%m%d_T%H%M%S')
    kwik_base = os.path.splitext(kwik_fn)[0]
    kwx_fn = kwik_base + '.kwx'


    with tb.open_file(kwik_fn, 'r') as kwik:

        try:
            o_dgrp = destination.get_node('/clusters')
            o_dgrp._f_move(newparent='/clusters_archive',
                           newname='clusters_%s' % o_dgrp._f_getattr('node_created_time'),
                           createparents=True,
                           )
            destination.flush()
        except tb.NoSuchNodeError:
            pass

        dgrp = destination.create_group('/', 'clusters')
        dgrp._f_setattr('node_created_time', dtg)
        dgrp._f_setattr('kwik_mod_time', 0.)  # we'll update this later when we find that everything goes as we plan.
        destination.flush()

        nshanks = kwik.root.channel_groups._v_nchildren
        nrecordings = kwik.root.recordings._v_nchildren
        rec_offsets = []
        for i in xrange(nrecordings):
            rec_grp = kwik.get_node('/recordings/%i' % i)
            rec_offsets.append(rec_grp._f_getattr('start_sample'))

        cluster_list = []
        for shank in xrange(nshanks):
            sh_grp = kwik.get_node('/channel_groups/%i'%shank)
            cluster_groups = {}
            for grp in sh_grp.cluster_groups.main:
                assert isinstance(grp, tb.Group)
                cluster_groups[grp._v_name] = grp._f_getattr('name')
            spike_times = sh_grp.spikes.time_samples.read()
            spike_recordings = sh_grp.spikes.recording.read()
            spike_clusters = sh_grp.spikes.clusters.main.read()

            # load spike waveforms if the kwx is present.
            if os.path.exists(kwx_fn):
                with tb.open_file(kwx_fn, 'r') as kwx:
                    kwx_grp = kwx.get_node('/channel_groups/%i' % shank)
                    sh_waveforms = kwx_grp.waveforms_filtered.read()
            else:
                sh_waveforms=None

            # correct for offset based on run starts. It is silly that kwik format doesn't just report spike
            # times relative to record 0, but this corrects that.
            for i in xrange(spike_times.size):
                rec = spike_recordings[i]  # the recording in which the spike occured specified by kwik
                spike_times[i] += rec_offsets[rec]

            for cluster in sh_grp.clusters.main:
                assert isinstance(cluster, tb.Group)
                cname = cluster._v_name
                cnum = int(cname)
                cgrp_num = cluster._f_getattr('cluster_group')
                try:
                    cgrp_name = cluster_groups[str(int(cgrp_num))]
                except TypeError as e:
                    print cgrp_num
                    print cluster
                    raise e
                if cgrp_name.lower() != 'noise':  # discard noise clusters.
                    # get the spike times:
                    clu_mask = spike_clusters == cnum
                    clu_spikes = spike_times[clu_mask]
                    if clu_spikes.size:  # if this is empty, don't add the cluster to the file.
                        carray = destination.create_carray(dgrp,
                                                           u's{0:02d}_c{1:04d}'.format(shank, cnum),
                                                           obj=clu_spikes)
                        carray.set_attr('cluster_group_num', cgrp_num)
                        carray.set_attr('cluster_group_name', cgrp_name)
                        carray.set_attr('shank', shank)
                        carray.flush()
                        cluster_list.append(carray)

                        # add mean waveform to the cluster attributes:
                        if sh_waveforms is not None:
                            clu_waveforms = sh_waveforms[clu_mask, :, :]
                            clu_mean_waveforms = np.mean(clu_waveforms, axis=0)
                            assert isinstance(carray, tb.CArray)
                            carray.set_attr('mean_waveforms', clu_mean_waveforms)

        for cluster in cluster_list:
            assert isinstance(cluster, tb.CArray)
            group_name = cluster._v_attrs['cluster_group_name']
            destination.move_node(cluster, '/clusters/{0}'.format(group_name), createparents=True)
        dgrp._f_setattr('kwik_mod_time', os.path.getmtime(kwik_fn))  # this will be updated to show that everything is ok.
    logging.info('cluster addition complete.')

    return