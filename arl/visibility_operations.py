# Tim Cornwell <realtimcornwell@gmail.com>
#
# Visibility data structure: a Table with columns ['uvw', 'time', 'antenna1', 'antenna2', 'vis', 'weight']
# and an attached attribute which is the frequency of each channel

import profile
import copy

from astropy import constants as const
from astropy.coordinates import SkyCoord, CartesianRepresentation
from astropy.table import Table, vstack

from crocodile.simulate import *

from arl.data_models import *
from arl.parameters import *

import logging

log = logging.getLogger("arl.visibility_operations")


def filter_gaintable(fg: GainTable, params={}):
    """Filter a Gaintable

    :param fg:
    :type GainTable:
    :returns: GainTable
    """
    # TODO: implement

    log.error("filter_gaintable: not yet implemented")
    return fg


def create_gaintable_from_array(gain: numpy.array, time: numpy.array, antenna: numpy.array, weight: numpy.array,
                                frequency: numpy.array, copy=False, meta=None, params={}):
    """ Create a gaintable from arrays

    :param gain:
    :type GainTable:
    :param time:
    :type numpy.array:
    :param antenna:
    :type numpy.array:
    :param weight:
    :type numpy.array:
    :param frequency:
    :type numpy.array:
    :param copy:
    :type bool:
    :param meta:
    :type dict:
    :param params: Dictionary containing parameters
    :returns: Gaintable
    """
    log_parameters(params)
    if meta is None:
        meta = {}
    nrows = time.shape[0]
    assert len(frequency) == gain.shape[1], "Discrepancy in frequency channels"
    assert len(antenna) == nrows, "Discrepancy in number of antenna rows"
    assert gain.shape[0] == nrows, "Discrepancy in number of gain rows"
    assert weight.shape[0] == nrows, "Discrepancy in number of weight rows"
    fg = GainTable()
    
    fg.data = Table(data=[gain, time, antenna, weight], names=['gain', 'time', 'antenna', 'weight'], copy=copy,
                    meta=meta)
    fg.frequency = frequency
    return fg


def interpolate_gaintable(gt: GainTable, params={}):
    """ Interpolate a GainTable to new sampling

    :param gt: GainTable
    :type GainTable:
    :param params: Dictionary containing parameters
    :returns: Gaintable
    """
    log_parameters(params)
    # TODO: implement

    log.error('"interpolate_gaintable: not yet implemented')
    return GainTable()


def combine_visibility(vis1: Visibility, vis2: Visibility, w1: float = 1.0, w2: float = 1.0, params={}) -> Visibility:
    """ Linear combination of two visibility sets

    :param vis1: Visibility set 1
    :type Visibility: Visibility to be processed
    :param vis2: Visibility set 2
    :type Visibility: Visibility to be processed
    :param w1: Weight of visibility set 1
    :type float:
    :param w2: Weight of visibility set 2
    :type float:
    :param params: Dictionary containing parameters
    :returns: Visibility
    """
    log_parameters(params)
    assert len(vis1.frequency) == len(vis2.frequency), "Visibility: frequencies should be the same"
    assert numpy.max(numpy.abs(vis1.frequency - vis2.frequency)) < 1.0, "Visibility: frequencies should be the same"
    assert len(vis1.data['vis']) == len(vis2.data['vis']), 'Length of output data table wrong'

    log.debug("visibility.combine: combining tables with %d rows" % (len(vis1.data)))
    log.debug("visibility.combine: weights %f, %f" % (w1, w2))

    v,w = combine_vis([vis1.vis, vis2.vis], [w1 * vis1.weight, w2 * vis2.weight])
    vis = Visibility(vis1, vis=v, weight=w)

    log.debug(u"combine_visibility: Created table with {0:d} rows".format(len(vis.data)))
    assert len(vis.data['vis']) == len(vis1.data['vis']), 'Length of output data table wrong'
    return vis


def concatenate_visibility(vis1: Visibility, vis2: Visibility, params={}) -> \
        Visibility:
    """ Concatentate the data sets in time, optionally phase rotating the second to the phasecenter of the first

    :param vis1:
    :type Visibility: Visibility to be processed
    :param vis2:
    :type Visibility: Visibility to be processed
    :param params: Dictionary containing parameters
    :returns: Visibility
    """
    log_parameters(params)
    assert len(vis1.frequency) == len(vis2.frequency), "Visibility: frequencies should be the same"
    assert numpy.max(numpy.abs(vis1.frequency - vis2.frequency)) < 1.0, "Visibility: frequencies should be the same"
    log.debug(
        "visibility.concatenate: combining two tables with %d rows and %d rows" % (len(vis1.data), len(vis2.data)))
    fvis2rot = phaserotate_visibility(vis2, vis1.phasecentre)
    vis = Visibility(vis1, data=vstack([vis1.data, fvis2rot.data], join_type='exact'))
    log.debug(u"concatenate_visibility: Created table with {0:d} rows".format(len(vis.data)))
    assert (len(vis.data) == (len(vis1.data) + len(vis2.data))), 'Length of output data table wrong'
    return vis


def flag_visibility(vis: Visibility, gt: GainTable = None, params={}) -> Visibility:
    """ Flags a visibility set, optionally using GainTable

    :param vis:
    :type Visibility: Visibility to be processed
    :param gt: GainTable
    :type GainTable:
    :param params: Dictionary containing parameters
    :returns: Visibility
    """
    # TODO: implement

    log_parameters(params)
    log.error("flag_visibility: not yet implemented")
    return vis


def filter_visibility(vis: Visibility, params={}) -> Visibility:
    """ Filter a visibility set

    :param vis:
    :type Visibility: Visibility to be processed
    :param params: Dictionary containing parameters
    :returns: Visibility
    """
    # TODO: implement

    log_parameters(params)
    log.error("filter_visibility: not yet implemented")
    return vis


def create_visibility(config: Configuration, times: numpy.array, freq: numpy.array, weight: float,
                      phasecentre: SkyCoord, meta: dict = None, params={}) -> Visibility:
    """ Create a Visibility from Configuration, hour angles, and direction of source

    :param config: Configuration of antennas
    :type Configuration:
    :param times: hour angles in radians
    :type numpy.array:
    :param freq: frequencies (Hz] Shape [nchan, npol]
    :type numpy.array:
    :param weight: weight of a single sample
    :type float:
    :param phasecentre: phasecentre of observation
    :type SkyCoord:
    :param meta:
    :type dict:
    :returns: Visibility
    """
    log_parameters(params)
    assert phasecentre is not None, "Must specify phase centre"
    nch = len(freq)
    npol = get_parameter(params, "npol", 4)
    ants_xyz = config.data['xyz']
    nants = len(config.data['names'])
    nbaselines = int(nants * (nants - 1) / 2)
    ntimes = len(times)
    nrows = nbaselines * ntimes
    row = 0
    rvis = numpy.zeros([nrows, nch, npol], dtype='complex')
    rweight = weight * numpy.ones([nrows, nch, npol])
    rtimes = numpy.zeros([nrows])
    rantenna1 = numpy.zeros([nrows], dtype='int')
    rantenna2 = numpy.zeros([nrows], dtype='int')
    for ha in times:
        rtimes[row:row + nbaselines] = ha * 43200.0 / numpy.pi
        for a1 in range(nants):
            for a2 in range(a1 + 1, nants):
                rantenna1[row] = a1
                rantenna2[row] = a2
                row += 1
    ruvw = xyz_to_baselines(ants_xyz, times, phasecentre.dec)
    log.debug(u"create_visibility: Created {0:d} rows".format(nrows))
    vis = Visibility()
    vis.data = Table(data=[ruvw, rtimes, rantenna1, rantenna2, rvis, rweight],
                     names=['uvw', 'time', 'antenna1', 'antenna2', 'vis', 'weight'], meta=meta)
    vis.frequency = freq
    vis.phasecentre = phasecentre
    vis.configuration = config
    return vis


def phaserotate_visibility(vis: Visibility, newphasecentre: SkyCoord, params={}) -> Visibility:
    """
    Phase rotate from the current phase centre to a new phase centre

    :param vis: Visibility to be rotated
    :type Visibility: Visibility to be processed
    :returns: Visibility
    """
    log_parameters(params)
    l, m, n = skycoord_to_lmn(newphasecentre, vis.phasecentre)
    log.debug('phaserotate_visibility: Relative cartesian representation of direction = (%f, %f, '
              '%f)' % (l, m, n))
    
    # Copy object and make a new table
    vis = copy.copy(vis)
    vis.data = vis.data.copy()
    
    # No significant change?
    if numpy.abs(l) > 1e-15 or numpy.abs(m) > 1e-15:
        log.debug('phaserotate: Phase rotation from %s to %s' % (vis.phasecentre, newphasecentre))
        
        # We are going to update in-place, so make a copy
        amp = vis.vis.copy()
        for channel in range(vis.nchan):
            uvw = vis.uvw_lambda(channel)
            phasor = simulate_point(uvw, l, m)
            for pol in range(vis.npol):
                log.debug('phaserotate: Phaserotating visibility for channel %d, polarisation %d' %
                          (channel, pol))
                amp[:, channel, pol] /= phasor
        
        # To rotate UVW, rotate into the global XYZ coordinate system and back
        xyz = uvw_to_xyz(vis.uvw, ha=-vis.phasecentre.ra, dec=vis.phasecentre.dec)
        uvw = xyz_to_uvw(xyz, ha=-newphasecentre.ra, dec=newphasecentre.dec)
    
        return Visibility(vis, vis=amp, uvw=uvw)

    return vis


def sum_visibility(vis: Visibility, direction: SkyCoord, params={}) -> numpy.array:
    """ Direct Fourier summation in a given direction

    :param vis: Visibility to be summed
    :type Visibility: Visibility to be processed
    :param direction: Direction of summation
    :type SkyCoord:
    :returns: flux[nch,npol], weight[nch,pol]
    """
    log_parameters(params)
    l,m,n = skycoord_to_lmn(direction, vis.phasecentre)
    log.debug('sum_visibility: Cartesian representation of direction = (%f, %f, %f)' % (
        l,m,n))
    flux = numpy.zeros([vis.nchan, vis.npol])
    weight = numpy.zeros([vis.nchan, vis.npol])
    for channel in range(vis.nchan):
        uvw = vis.uvw_lambda(channel)
        phasor = numpy.conj(simulate_point(uvw, l,m))
        for pol in range(vis.npol):
            log.debug('sum_visibility: Summing visibility for channel %d, polarisation %d' % (
                channel, pol))
            ws = vis.weight[:, channel, pol]
            wvis = ws * vis.vis[:, channel, pol]
            flux[channel, pol] += numpy.sum(numpy.real(wvis * phasor))
            weight[channel, pol] += numpy.sum(ws)
    flux[weight > 0.0] = flux[weight > 0.0] / weight[weight > 0.0]
    flux[weight <= 0.0] = 0.0
    return flux, weight


def coalesce_visibility(vis: Visibility, params={}) -> Visibility:
    """ Coalesce visibilities in time and frequency according to baseline length
    
    Creates new Visibility by averaging in time and frequency
    
    :param vis: Visibility to be coalesced
    :type Visibility:
    :returns: Visibility after coalescing
    """
    # TODO: implement

    log_parameters(params)
    log.error("coalesce_visibility: not yet implemented")
    return vis


def de_coalesce_visibility(vis: Visibility, vistemplate: Visibility, params={}) -> Visibility:
    """ De-coalesce visibility in time and frequency i.e. replicate to template Visibility
    
    This is the opposite of coalescing - the Visibility is expanded into sampling independent
    of baseline length.
    
    :param vis: Visibility to be de-coalesced
    :type Visibility: Visibility
    :param vistemplate: template Visibility
    :type Visibility: Visibility
    :returns: Visibility after de-coalescing
    """
    # TODO: implement

    log_parameters(params)
    log.error("de_coalesce_visibility: not yet implemented")
    return vis


def aq_visibility(vis, params={}):
    """Assess the quality of Visibility

    :param vis: Visibility to be assessed
    :type Visibility:
    :returns: AQ
    """
    context = get_parameter(params, 'context', None)
    log_parameters(params)
    avis = numpy.abs(vis.vis)
    data = {'maxabs': numpy.max(avis),
            'minabs': numpy.min(avis),
            'rms': numpy.std(avis),
            'medianabs': numpy.median(avis)}
    qa = QA(origin=None,
            data=data,
            context=get_parameter(params, 'context', None))
    return qa
