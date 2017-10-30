#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ROTATIONAL SEISMOLOGY ROUTINES. Following the theory and assuming a
transversely polarized plane wave, this script compares the transverse 
acceleration and vertical rotation rate of an event through:

1) direct waveform comparison in different time windows (P-coda, S-waves and
surface waves) using broadband data from WET (Wettzell) station for 
acceleration and from Wettzell ringlaser (RLAS station) for the vertical
rotation rate.

2) zero-lag cross-correlation analysis of the waveforms using the theoretical
backazimuth of the event. Aside from the correlation-coefficients, it returns
the local horizontal phase velocity for the time windows with a correlation
factor larger than 0.75. It also preforms correlation analysis for different
azimuths and estimates the backazimuth from the one that has the highest
correlation coefficients.

Additionally, the routine generates:

+ a QuakeML (.xml) file that stores data for each event, and contains an extra 
tag with rotational parameters, to be used in the JANE format database

+ a human readable (.json) ordered dictionary text file that contains both
event and station information, as well as output data results from processing


INFORMATION: This script can be used to generate figures for the database.

+ QuakeML files can be read in directly, for example to pull events from other 
catalogs (ISC, ...) by setting: --mode qmlfile

+ Events can be chosen from the IRIS FDSN catalog, which is usually faster
and more abundant especially for older events. Set: --mode fdsn

+ Events are bandstoppped for the secondary microseism (5-12s) if they are 
non-local.

+ maximum correlation coefficients for the estimated BAz are added 
in the subplot 3 on page 3.

+ only 0.5 hour recordings shown for local and close events!

+ P-coda windows (sec_p) now shorter for local and close events (2s).

+ 28.09.17 - To-change log: 
    -Get_corrcoeffs: clean up the cross correlation for loops
    -phas_vel: numpy mean of empty slice throwing runtime warning
    -filter_and_rotate: can compress the different filter bands into a few lines
    -resample: no else finish - may throw error?
    -plot_wa...: figure out what to do with the flip call
        -minimize misfit or waterlevel method at taper
    -phase velocities:
        threshold value of corrcoef as a parameter
        phase velocity estimation by misfit minimization

"""
import os
import sys
import json
import glob
import obspy
import shutil
import argparse
import datetime
import warnings
from collections import OrderedDict
from urllib.request import urlopen
from xml.dom.minidom import parseString

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from obspy.core import read
import matplotlib.pylab as plt
from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy import read_events, Catalog
from mpl_toolkits.basemap import Basemap
# from obspy.imaging.beachball import Beach
from obspy.signal.rotate import rotate_ne_rt
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.util.attribdict import AttribDict
from obspy.clients.fdsn import Client as fdsnClient
from obspy.signal.cross_correlation import correlate
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees

# warnings.filterwarnings(
#     action='once', category=np.VisibleDeprecationWarning,
#     message='using a non-integer number instead of an integer will result in '
#             'an error in the future')

# if matplotlib.__version__ < '1.0':  # Matplotlib 1.0 or newer is necessary
#     raise ValueError('I need Matplotlib version 1.0 or newer.')

class RotationalProcessingException(Exception):
    """
    Exception for when no data can be found for an event
    """
    pass

def download_data(origin_time, instrument_id, source):
    """
    It downloads the data from seismic stations for the desired event(s).
    Inputs are the origin time (UTC), network, station, location and channel
    of the event. Returns a stream object fetched from Arclink. If Arclink
    does not work data is alternatively fetched from Seishub.

    :type origin_time: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param origin_time: origin time of the event.
    :type net: str
    :param net: Network code, e.g. ``'BW'``.
    :type sta: str
    :param sta: Station code, e.g. ``'WET'``.
    :type loc: str
    :param loc: Location code, e.g. ``'01'``. Location code may
        contain wild cards.
    :type chan: str
    :param chan: Channel code, e.g. ``'EHE'``. Channel code may
        contain wild cards.

    :type st: Stream object :class: `~obspy.core.stream.Stream`
    :return st: fetched data stream
    """
    # arclink is deprecated, but call is kept for posterity? rip
    # try:
    #     c = arclinkClient(user='test@obspy.org')
    #     st = c.get_waveforms(network=net, station=sta, location='', 
    #                          channel=chan,
    #                          starttime=origin_time-190,
    #                          endtime=origin_time+3*3600+10)
    
    # check paths to see if running on FFB, LMU or neither
    st = None
    dataDir_get = '/bay200/mseed_online/archive/' #FFB
    if not os.path.exists(dataDir_get):
        dataDir_get = '/import/netapp-m-02-bay200/mseed_online/archive/'#LMU            
    if not os.path.exists(dataDir_get):
        dataDir_get = None
    
    net, sta, loc, cha = instrument_id.split('.')
    
    # if data path exists, read in data from file
    if dataDir_get:
        print("Fetching {} data from file".format(net))
        fileName = '.'.join(instrument_id,'D',origin_time.strftime('%Y.%j'))
        filePath = os.path.join(dataDir_get, origin_time.strftime('%Y'),
                                net, sta, cha + '.D', fileName)
        
        origin_time2 = origin_time + 86400
        fileName2 = '.'.join(instrument_id,'D',origin_time2.strftime('%Y.%j'))
        filePath2 = os.path.join(dataDir_get, origin_time2.strftime('%Y'),
                                net, sta, chan + '.D', fileName2)

        if os.path.isfile(filePath):
            data_source = 'Archive'
            if origin_time.hour > 21:
                st = Stream()
                st.extend(read(filePath, starttime = origin_time - 180,
                      endtime = origin_time + 3 * 3600))
                st.extend(read(filePath2, 
                      starttime = UTCDateTime(o_time2.year, o_time2.month, 
                                                            o_time2.day, 0, 0),
                      endtime = origin_time + 3 * 3600))
                st.merge(method=-1)
            else:
                st = read(filePath, starttime = origin_time - 180,
                                    endtime = origin_time + 3 * 3600)    
        else:
            print("\tFile not found: \n\t %s \n" % filePath)    
    
    # if data/path does not exist, try querying FDSN webservices
    elif (not dataDir_get) or (not st):
        for S in source:
            try:
                print("Fetching {} data from FDSN ({})".format(net,S))
                c = fdsnClient(S)
                st = c.get_waveforms(network=net, station=sta, location=loc, 
                                    channel=cha, starttime=origin_time-190,
                                    endtime=origin_time+3*3600+10)
                break
            except:
                print("\tFailed")
                pass
        data_source = S 
    
    if not st:
        raise RotationalProcessingException(
                                        "Data not available for this event")

    st.trim(starttime=origin_time-180, endtime=origin_time+3*3600)
    print("\tDownload of {!s} {!s} data successful".format(
              st[0].stats.station, st[0].stats.channel))
    
    return st, data_source


def event_info_data(event, station, mode, polarity, instrument):

    """
    Extracts information from the event and generates variables containing
    the event latitude, longitude, depth, and origin time.
    Ringlaser (RLAS) and broadband signals (WET) are received from the
    download_data function.
    The great circle distance (in m and °) between event location and station
    in Wetzell, as well as the theoretical backazimuth are computed.

    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :type mode: str
    :param mode: Defines where data fetched from
    :type polarity: str
    :param polarity: 'normal' or 'reverse' choice for rotation polarity
    :type instrument: str
    :param instrument: 'wet' or 'wetr' choice for translation data
    :rtype event_lat: float
    :return event_lat: Latitude of the event in degrees.
    :rtype event_lon: float
    :return event_lon: Longitude of the event in degrees.
    :rtype depth: float
    :return depth: Hypocenter depth in km
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :return startev: Origin time of the event.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Three component broadband station signal.
    :rtype baz: tuple
    :return baz: [0] great circle distance in m, 
                 [1] theoretical azimuth,
                 [2] theoretical backazimuth.
    """
    origin = event.preferred_origin() or event.origins[0]
    startev = origin.time

    if station == 'RLAS':
        station_lat = 49.144001
        station_lon = 12.8782

        source = ['http://eida.bgr.de', 
                  'http://erde.geophysik.uni-muenchen.de']

        # ringlaser signal, source LMU first
        rotation_id = 'BW.RLAS..BJZ'
        if origin.time < UTCDateTime(2010, 4, 16):
            rotation_id = 'BW.RLAS..BAZ' 

        rt,srcRT = download_data(startev, rotation_id, source[::-1])
        if polarity == 'reverse':
            rt[0].data *= -1

        # create dictionary for sources
        data_sources = {'BJZ':srcRT}

        # broadband station signal, assume all translation same source
        ac = Stream()
        for channels in ['BHN','BHE','BHZ']:
            if instrument == 'STS2':
                translation_id = 'GR.WET..{}'.format(channels)
            elif instrument == 'LENNARTZ':
                translation_id = 'BW.WETR..{}'.format(channels)

            tr,srcTR = download_data(startev, translation_id, source)
            data_sources[channels] = srcTR
            ac += tr

    elif station == 'ROMY':
        station_lat = 48.162941
        station_lon = 11.275476

        source = ['http://eida.bgr.de', 
                  'http://erde.geophysik.uni-muenchen.de']

        # ringlaser signal, source LMU first
        rotation_id = 'BW.ROMY..BJZ'
        rt,srcRT = download_data(startev, rotation_id, source[::-1])
        if polarity.lower() == 'reverse':
            rt[0].data *= -1

        # create dictionary for sources
        data_sources = {'BJZ':srcRT}

        # broadband station signal, assume all translation has source
        ac = Stream()
        for channels in ['BHN','BHE','BHZ']:
            translation_id = 'GR.FUR..{}'.format(channels)
            tr,srcTR = download_data(startev, translation_id, source)
            data_sources[channels] = srcTR
            ac += tr
        
    # coordinates, theoretical event backazimuth and distance
    event_lat = origin.latitude
    event_lon = origin.longitude
    depth = origin.depth * 0.001  # Depth in km
    dist_baz = gps2dist_azimuth(lat1 = event_lat, lon1 =  event_lon, 
                           lat2 = station_lat, lon2 = station_lon)

    for ca in [ac[0], ac[1], ac[2], rt[0]]:
        ca.stats.coordinates = AttribDict()
        ca.stats.coordinates['longitude'] = station_lat
        ca.stats.coordinates['latitude'] = station_lon
        ca.stats['back_azimuth'] = dist_baz[2]
        ca.stats['starttime'] = startev - 180
        ca.stats['sampling_rate'] = 20.
    
    return event_lat, event_lon, depth, startev, rt, ac, dist_baz, data_sources


def is_local(distance):

    """
    Checks whether the event is close (< 333.33 km), local (< 1111.1 km) or
    non-local.

    :type distance: float
    :param distance: event station distance in km
    :rtype: str
    :return: Self-explaining string for event distance.
    """
    if 0.001 * distance / 111.11 < 10.0:
        if 0.001 * distance / 111.11 < 3.0:
            is_local = 'close'
        else:
            is_local = 'local'
    else:
        is_local = 'non-local'

    return is_local


def get_moment_tensor_magnitude(link):

    """
    Extracts the moment tensor and magnitude for an event if an IRIS-xml file
    is given by a url (link).

    :type link: str
    :param link: Link to the IRIS-xml, where event and moment tensor data are
        fetched.
    :rtype MomentTensor: list of floats
    :return MomentTensor: List of the six independent components of the moment
        tensor.
    :rtype Magnitude: float
    :return Magnitude: Moment magnitude (Mw) of the event.
    """
    file = urlopen(link)
    data = file.read()
    file.close()

    dom = parseString(data)
    xmlTag = dom.getElementsByTagName('value')
    xmlTag2 = dom.getElementsByTagName('text')
    
    # 19th value in xml-file
    Magnitude = float(xmlTag[19].firstChild.nodeValue)
    Region = str(xmlTag2[0].firstChild.nodeValue)
    MomentTensor = []

    for i in range(1, 7):
        value = float(xmlTag[i].firstChild.nodeValue)
        MomentTensor.append(value)

    return MomentTensor, Magnitude, Region


def resample(is_local, rt, ac):

    """
    Resamples signal accordingly with sampling rates and cut-off frequencies
    dependent on the location of the event (5 sec and 2Hz for local events,
    60 sec and 1 Hz for non-local events).

    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :type baz: tuple
    :param baz: [0] Great circle distance in m, 
                [1] azimuth A->B in degrees,
                [2] backazimuth B->A in degrees.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Decimated rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Decimated three component broadband station signal.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: (Decimated) copy of rotational signal from ringlaser
        for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: (Decimated) copy of the three component broadband
        station signal for p-coda calculation.
    :rtype sec: int
    :return sec: Time window length.    
    :rtype cutoff: float
    :return cutoff: Cut-off frequency for the lowpass filter.
    :rtype cutoff_pc: float
    :return cutoff_pc: Cut-off frequency for the highpass filter in P-coda.
    """

    cutoff_pc = 0.5  # Cut-off frequency for the highpass filter in the P-coda
    if is_local == 'non-local':
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt_pcoda.decimate(factor=2)
        ac_pcoda.decimate(factor=2)
        rt.decimate(factor=4)
        ac.decimate(factor=4)
        sec = 120 # Length of time window in seconds
        sec_p = 5 # Pcoda time window in seconds
        cutoff = 1.0 # Cut-off freq for the lowpass filter
    elif is_local == 'local':
        for trr in (rt + ac):
            trr.data = trr.data[0: int(1800 * rt[0].stats.sampling_rate)]
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 5 
        sec_p = 2
        cutoff = 2.0  
    elif is_local == 'close':
        for trr in (rt + ac):
            trr.data = trr.data[0: int(1800 * rt[0].stats.sampling_rate)]
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 3
        sec_p = 2
        cutoff = 4.0  

    return rt, ac, rt_pcoda, ac_pcoda, sec, sec_p, cutoff, cutoff_pc


def remove_instr_resp(rt, ac, rt_pcoda, ac_pcoda, station, startev):

    """
    This function removes the instrument response from the original signal
    and checks if starttime and endtime match for both instruments. 
    * Note for Poles and Zeros:
        1 zero acceleration/ 2 zeros velocity/ 3 zeros displacement for STS2 


    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rotational signal from ringlaser
        for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of the three component broadband
        station signal for p-coda calculation.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param startev: Origin time of the event.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Detrended and trimmed rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Detrended and trimmed three component broadband station signal.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: Detrended and trimmed copy of rotational signal from
        ringlaser for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: Detrended and trimmed copy of the three component
        broadband station signal for p-coda calculation.
    """

    # poles and zeros dictionaries for units of nanometers/s^2
    paz_sts2 = {'poles': [(-0.0367429 + 0.036754j),
                          (-0.0367429 - 0.036754j)],
                'sensitivity': 0.944019640, 
                'zeros': [0j], 
                'gain': 1.0}
    paz_lennartz = {'poles': [(-0.22 + 0.235j),
                                (-0.22 - 0.235j), 
                                (-0.23 + 0.0j)],
                    'sensitivity': 1, 
                    'zeros': [(0+0j),(0+0j)], 
                    'gain': 1.0}

    # standard preprocessing 
    rt.detrend(type='linear')
    ac.detrend(type='linear')
    ac.taper(max_percentage=0.05)
    rt.taper(max_percentage=0.05)
    ac_pcoda.detrend(type='linear')
    rt_pcoda.detrend(type='linear')

    # remove instrument response
    if station == 'RLAS':
        rt[0].data = rt[0].data * 1/(6.3191e3)  # rotation rate: nrad/s
        rt_pcoda[0].data = rt_pcoda[0].data * 1/(6.3191e3)

        if instrument == 'STS2':
            ac.simulate(paz_remove=paz_sts2, remove_sensitivity=True)  # nm/s^2
            ac_pcoda.simulate(paz_remove=paz_sts2, remove_sensitivity=True)

        elif instrument == 'LENNARTZ':
            ac.simulate(paz_remove=paz_lennartz, remove_sensitivity=True)
            ac_pcoda.simulate(paz_remove=paz_lennartz, remove_sensitivity=True)
        
            ac.filter('highpass', freq=0.04, zerophase=True, corners=3)
            ac_pcoda.filter('highpass', freq=0.04, zerophase=True, corners=3)         
     
           
    elif station == 'ROMY':
        rt[0].data = rt[0].data * 1/(1.01821e4) # rotation rate in nrad/s
        rt_pcoda[0].data = rt_pcoda[0].data * 1/(1.01821e4)  # nrad/s
 
        ac.simulate(paz_remove=paz_sts2, remove_sensitivity=True)  # nm/s^2
        ac_pcoda.simulate(paz_remove=paz_sts2, remove_sensitivity=True)


    elif station == 'PFO':
        rt[0].data = rt[0].data * 1. / 2.5284 * 1e-3  # rotation rate in nrad/s
        rt_pcoda[0].data = rt_pcoda[0].data * 1. / 2.5284 * 1e-3  # nrad/s
        ac.remove_response(output='ACC', pre_filt=(0.005, 0.006, 30., 35.))
        ac_pcoda.remove_response(output='VEL',
                                 pre_filt=(0.005, 0.006, 30., 35.))

        # to nm/s^2
        for traza in (ac + ac_pcoda):
            traza.data = 1e9 * traza.data

    else:
        sys.exit('Invalid station')

    # make sure start and endtimes match for both instruments, if not, trim
    startaim = max([tr.stats.starttime for tr in (ac + rt)])
    endtaim = min([tr.stats.endtime for tr in (ac + rt)])
    ac.trim(startaim, endtaim, nearest_sample=True)
    rt.trim(startaim, endtaim, nearest_sample=True)
    ac_pcoda.trim(startaim, endtaim, nearest_sample=True)
    rt_pcoda.trim(startaim, endtaim, nearest_sample=True)

    return rt, ac, rt_pcoda, ac_pcoda

def gaussian_filter(sigarray, delta, bandwidth, freq0):
    
    """
    Gaussian filter function. 
    Not used in script, left incase necessary later
    
    :type sigarray: np array
    :param sigarray: signal array (much faster if the length is a power of 2)
    :type delta: float
    :param delta: time sampling interval (seconds)
    :type bandwidth: float    
    :param bandwith: filter df (must be > 0)
    :type freq0: float
    :param freq0: center frequency (Hz)
    :rtype sigarray_filtered: np array
    :return sigarray_filtered: the guassian filtered array

    """
    # prepare the frequency domain filter
    n = len(sigarray)  # N
    freq = fftfreq(n, delta) # exact (!) frequency array
    
    # construct our gaussian according the constQ criterion of Archambeau et al.
    # do not forget negative frequencies
    beta = np.log(2.)/2.
    g = np.sqrt(beta/np.pi)*np.exp(-beta * 
                                    (np.abs(freq - freq0) / bandwidth) ** 2.) 


    # convolve your signal by the filter in frequency domain 
    sigarray_fourier = fft(sigarray) 
    sigarray_fourier_filtered = sigarray_fourier * g

    # back to time domain
    sigarray_filtered = np.real(ifft(sigarray_fourier_filtered))
    sigarray_filtered = highpass(sigarray_filtered, freq=0.0033, 
                                                df=5, corners=3, zerophase=True)
    sigarray_filtered = detrend(sigarray_filtered)

    return sigarray_filtered


def filter_and_rotate(ac, rt, baz, rt_pcoda, ac_pcoda, cutoff, cutoff_pc,
                      station, is_local):


    """
    Filters trace data using the cut-off frequencies and
    lowpass/ highpass and zerophase filters. Rotates the horizontal components
    of the signal to theoretical backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rotational signal from ringlaser
        for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of the three component broadband
        station signal p-coda calculation.
    :type cutoff: float
    :param cutoff: Cut-off frequency for the lowpass filter.
    :type cutoff_pc: float
    :param cutoff_pc: Cut-off frequency for the highpass filter in P-coda.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :rtype rotate: :class: `~obspy.core.stream.Stream`
    :return rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :rtype pcod_rotate: :class: `~obspy.core.stream.Stream`
    :return pcod_rotate: Stream object of the broadband station signal with
        rotated horizontal components for P-coda calculations.
    :rtype pcoda_rotate: :class: `~obspy.core.stream.Stream`
    :return pcoda_rotate: Stream object of the broadband station signal with
        rotated horizontal components for P-coda calculations.
    :rtype frtp: :class: `~obspy.core.stream.Stream`
    :return frtp: Filtered rt_pcoda.
    :rtype facp: :class: `~obspy.core.stream.Stream`
    :return facp: Filtered ac_pcoda.
    :rtype frotate: :class: `~obspy.core.stream.Stream`
    :return frotate: Filtered and rotated ac_pcoda.
    :rtype cop_rt: :class: `~obspy.core.stream.Stream`
    :return cop_rt: Highpass filtered rotational trace (rt).
    """

    # set the list of frequencies for to bandpass filter for
    freq_list = [0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0]
    number_of_bands = len(freq_list) - 1
    
    # copies for phase velocities and pcoda analysis
    ac_pcoda_coarse = ac.copy()
    rt_pcoda_coarse = rt.copy()
    rt_bands = [rt.copy() for _ in range(number_of_bands)]
    rotate_bands = [ac.copy() for _ in range(number_of_bands)]

    # filter streams high and low
    ac.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    rt.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    ac.filter('highpass', freq=0.005, corners=2, zerophase=True)
    rt.filter('highpass', freq=0.005, corners=2, zerophase=True)
    ac_pcoda_coarse.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    rt_pcoda_coarse.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)

    # filter out secondary microseisms (5-12s) for non-local events
    if is_local == "non-local":    
        ac.filter('bandstop', freqmin=1/12, freqmax=1/5, 
                                corners=4, zerophase=True)
        rt.filter('bandstop', freqmin=1/12, freqmax=1/5, 
                                corners=4, zerophase=True)

    # rotate translational signals to theoretical event backazimuth
    rotate = ac.copy()
    pcoda_rotate = ac_pcoda.copy()
    pcoda_rotate_coarse = ac_pcoda_coarse.copy()

    rotate.rotate(method = 'NE->RT')
    pcoda_rotate.rotate(method = 'NE->RT')
    pcoda_rotate_coarse.rotate(method = 'NE->RT')
    
    # for varying frequency bands, rotate to BAz, filter both RT and Rotate 
    for I in range(number_of_bands):
        rotate_bands[I].rotate(method='NE->RT')

        for band_select in [rt_bands[I],rotate_bands[I]]:
            band_select[0].filter(type = 'bandpass',
                                freqmin = freq_list[I],
                                freqmax = freq_list[I+1],
                                corners = 3,
                                zerophase = True)

    # Filter Pcoda for BAz analysis
    filt_rt_pcoda = rt_pcoda.copy()
    filt_ac_pcoda = ac_pcoda.copy()
    filt_rotate = filt_ac_pcoda.copy()
    filt_rt_pcoda.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    filt_ac_pcoda.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    filt_rotate.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    filt_rotate.rotate(method = 'NE->RT')


    return rt_bands, rotate_bands, rotate, rt_pcoda_coarse, pcoda_rotate_coarse, \
            filt_rt_pcoda, filt_ac_pcoda, filt_rotate 



def ps_arrival_times(distance, depth, init_sec):

    """
    Obtains the arrival times (in seconds after the start time of the fetched
    data) of the first P an S waves of the event. The inputs are the
    epicentral distance in degrees, the depth in km and the initial time in
    seconds (starttime_of_the_event - data_starttime)

    :type distance: float
    :param distance: Great circle distance between earthquake source and
        receiver station.
    :type depth: float
    :param depth: Hypocenter depth in km.
    :type init_sec: float
    :param init_sec: Initial time of the event in sec in the fetched data.
    :rtype arriv_p: float
    :return arriv_p: Arrival time of the first P-wave.
    :rtype arriv_s: float
    :return arriv_s: Arrival time of the first S-wave.
    """
    # use taup to get the theoretical arrival times for P & S
    TauPy_model = TauPyModel('iasp91')
    tt = TauPy_model.get_travel_times(
                                distance_in_degree=0.001 * distance / 111.11, 
                                source_depth_in_km=depth)
    
    times_p,times_s = [],[]
    # from all possible P arrivals, select the earliest 
    pwave_list = ['P','p','Pdiff','PKiKP','PKIKP','PP','Pb','Pn','Pg']
    for i2 in range(len(tt)):
        if tt.__getitem__(i2).__dict__['name'] in pwave_list:
            ptime = tt.__getitem__(i2).__dict__['time']
            times_p.append(ptime)


    # from all possible S arrivals, select the earliest
    swave_list = ['S','s','Sdiff','SKiKS','SKIKS','SS','Sb','Sn','Sg']
    for i3 in range(len(tt)):
        if tt.__getitem__(i3).__dict__['name'] in swave_list:
            stime = tt.__getitem__(i3).__dict__['time']
            times_s.append(stime)

    arriv_p = np.floor(init_sec + min(times_p))
    arriv_s = np.floor(init_sec + min(times_s))

    return arriv_p, arriv_s


def time_windows(DS, arriv_p, arriv_s, init_sec, is_local):
    """
    Determines time windows for arrivals and subplots for P-waves,
    S-waves, initial and latter surface waves. 
    Window lengths dependent on event distance.
    Keep all values integers, as they are used for slice indexing



    :type DS: float
    :param DS: event station distance in km
    :type arriv_p: float
    :param arriv_p: Arrival time of the first P-wave.
    :type arriv_s: float
    :param arriv_s: Arrival time of the first S-wave.
    :type init_sec: float
    :param init_sec: Initial time of the event in sec in the fetched data.
    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :rtype: int
    :return min_pw: Starttime for P-waves window.
    :return max_pw: Endtime for P-waves window.
    :return min_sw: Starttime for S-waves window.
    :return max_sw: Endtime for S-waves window.
    :return min_lwi: Starttime for initial surface-waves window.
    :return max_lwi: Endtime for initial surface-waves window.
    :return min_lwf: Starttime for latter surface-waves window.
    :return max_lwf: Endtime for latter surface-waves window.
    """

    if is_local == 'non-local':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + (arriv_s - arriv_p) // 4)
        min_sw = int(round(arriv_s - 0.001 * (arriv_s - arriv_p)))
        max_sw = int(arriv_s + 150)
        min_lwi = int(round(surf_tts(DS, init_sec) - 20))
        max_lwi = int(min_lwi + round((DS/1E6) * 50)) # 50 sec per 1000 km. 
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + round((DS/1E6) * 60)) # 60 sec per 1000 km.
    elif is_local == 'local':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + 20)
        min_sw = int(arriv_s - 5)
        max_sw = int(min_sw + 20)
        min_lwi = int(round(surf_tts(DS, init_sec) + 20))
        max_lwi = int(min_lwi + 50)
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + 80)
    elif is_local == 'close':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + 7)
        min_sw = int(arriv_s)
        max_sw = int(min_sw + 7)
        min_lwi = int(round(surf_tts(DS, init_sec) + 5))
        max_lwi = int(min_lwi + 12)
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + 80)

    # I wrote a dictionary but it's not so useful... :/

    # arrival_dic = {'pwave_start': min_pw,
    #             'pwave_end': max_pw,
    #             'swave_start': min_sw,
    #             'swave_end': max_sw,
    #             'initial_surface_start': min_lwi,
    #             'initial_surface_end': max_lwi,
    #             'latter_surface_start': min_lwf,
    #             'latter_surface_end': max_lwf}


    return min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf


def surf_tts(distance, start_time):

    """
    Uses arrival times for different epicentral distances based on the IASP91
    travel times model to estimate a curve of travel times for surface waves
    and get the arrival time of the surface waves of the event. Inputs are the
    epicentral distance in degrees and the event start time in seconds.

    :type distance: float
    :param distance: Epicentral distance in degrees between earthquake source
        and receiver station.
    :type start_time: float
    :param start_time: Starttime of the event in the fetched seismogram.
    :rtype arrival: float
    :return arrival: Arrival time of the surface waves of the event.
    """
    deltas = np.arange(0., 140., 5.)
    tts = 60. * np.array(
        [0., 2., 4., 6.2, 8.4, 11., 13., 15.2, 17.8, 19.4, 22., 24.1, 26.6,
         28.6, 30.8, 33., 35.6, 37.4, 39.8, 42., 44.2, 46.4, 48.8, 50.9, 53.6,
         55.2, 57.8, 60.])
    (mval, nval) = np.polyfit(deltas, tts, 1)

    # calculate surface wave travel times for degrees 1 to 180 ?
    surftts = mval * np.arange(0., 180.1, 0.01)
    difer = []
    for i4 in range(0, len(surftts)):
        dife_r = abs(0.001 * distance / 111.11 - np.arange(0., 180.1, 0.01)
                     [i4])
        difer.append(dife_r)
    # love wave arrival: event time + surftts for closest degree??
    # (smallest difference between distance for surftts and actual distance of
    #  event)

    arriv_lov = np.floor(start_time + surftts[np.asarray(difer).argmin()])
    diferans = []
    for i1 in range(len(deltas)):
        dif2 = abs(np.arange(0., 180.1, 0.01)[np.asarray(difer).argmin()] -
                   deltas[i1])
        diferans.append(dif2)
    # arrival = love wave arrival - p arrival?
    peq = surftts[np.asarray(difer).argmin()] - \
        tts[np.asarray(diferans).argmin()]
    arrival = arriv_lov + peq

    return arrival


def get_corrcoefs(streamA, streamB, sec, station):

    """
    Calculates the zero-lag correlation coefficients between two streams
    *StreamA and StreamB data lengths need to be the same.

    :type rt: :class: `~obspy.core.trace.Trace`
    :param rt: Trace of the rotational data from ringlaser.
    :type rodat: numpy.ndarray
    :param rodat: Rotational data ...
    :type acstr: :class: `~obspy.core.stream.Stream`
    :param acstr: Three component broadband station signal.
    :type rotate_array: numpy.ndarray
    :param rotate_array:
    :type sec: int
    :param sec: Time window length.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :rtype corrcoefs: numpy.ndarray
    :return corrcoefs: Correlation coefficients.
    :rtype thres: numpy.ndarray
    :return thres: Array for plotting dashed line of 75'%' correlation.
    """

    # time window in samples
    strA_TW = int(streamA[0].stats.sampling_rate * sec)
    strB_TW = int(streamB[0].stats.sampling_rate * sec)

    corrcoefs = []
    for i in range(0, len(streamA[0].data) // strA_TW):
        coeffs = correlate(a = streamA[0].data[i*strA_TW:(i+1)*strA_TW],
                           b = streamB[0].data[i*strB_TW:(i+1)*strB_TW], 
                           shift=0)

        corrcoefs.append(coeffs[0])

    corrcoefs = np.asarray(corrcoefs)
    thres = 0.75 * np.ones(len(corrcoefs) + 1)

    return corrcoefs, thres


def baz_analysis(rt, ac, sec, station):

    """
    Backazimuth analysis: Computes the correlation coefficients for
    the backazimuth and backazimuth values.
    For each time window, loop through all backazimuth values, calculate 
    correlations for each backazimuth, iterate over entire trace time
    *streamA and streamB need to be the same length

    :type rt: :class: `~obspy.core.trace.Trace`
    :param rt: Trace of the rotational data from ringlaser.
    :type rodat: numpy.ndarray
    :param rodat: Rotational data ...
    :type acstr: :class: `~obspy.core.stream.Stream`
    :param acstr: Three component broadband station signal.
    :type sec: int
    :param sec: Time window length.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type ind: int
    :param ind: Index for stream data selection.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :rtype corrbaz: numpy.ndarray
    :return corrbaz: Correlation coefficients for each backazimuth.
    :rtype maxcorr: numpy.ndarray
    :return maxcorr: Backazimuth values for maximum correlation for each time
        window.
    :rtype backas: numpy.ndarray
    :return backas: Vector containing backazimuths (step: 10°).
    """

    # time window in samples
    rt_TW = int(rt[0].stats.sampling_rate * sec)
    ac_TW = int(ac[0].stats.sampling_rate * sec)

    # create a list of backazimuths to iterate over
    step = 10
    corr_length = len(rt[0].data) // rt_TW
    backas = np.linspace(0, 360 - step, 360 / step)
    
    # iterate over BAz, rotate, correlate trace
    corrbaz_list = []
    for BAZ in backas:
        for j in range(0, corr_length):
            rotate_tmp = ac.copy()
            rotate_tmp.rotate(method='NE->RT', back_azimuth = BAZ)
            trv_acc = rotate_tmp.select(component='T')

            corrbaz = correlate(a = rt[0].data[j*rt_TW:(j+1)*rt_TW],
                                b = trv_acc[0].data[j*ac_TW:(j+1)*ac_TW],
                                shift = 0)
            corrbaz_list.append(corrbaz[0])

    corrbaz_list = np.asarray(corrbaz_list)
    corrbaz_list = corrbaz_list.reshape(len(backas), corr_length)

    # find maximum correlations
    maxcorr_list = []
    for k in range(0, corr_length):
        maxcorr_list.append(backas[corrbaz_list[:,k].argmax()])

    maxcorr_list = np.asarray(maxcorr_list)

    # array containing max corr coef for each window
    coefs_list = []    
    for l in range(0, corr_length):
        coefs_list.append(np.max(corrbaz_list[:,l]))

    return corrbaz_list, maxcorr_list, backas, coefs_list


def estimate_baz(rt, trv_acc, start_sw, end_sw, station):
    """
    Calculates the sum of all correlation coefficients above a certain
    threshold (0.9) within S-waves and surface waves for each backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type min_sw: float
    :param min_sw: Starttime for S-waves window.
    :type max_lwf: float
    :param max_lwf: Endtime for latter surface waves window.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :rtype corrsum: list of floats
    :return corrsum: Sum of all correlation coefficients above a certain
        threshold (0.9) within S-waves and surface waves for each backazimuth.
    :rtype baz_list: numpy.ndarray
    :return baz_list: Vector containing backazimuths (step: 1°).
    """

    # set integer sampling rates
    rt_SR = int(rt[0].stats.sampling_rate)
    trv_SR = int(trv_acc[0].stats.sampling_rate)
    
    # time window in samples
    sec_internal = 30
    rt_TW = sec_internal * rt_SR
    trv_TW = sec_internal * trv_SR

    # sample number of surface wave start/end
    start_sample = start_sw * rt_SR
    end_sample = end_sw * trv_SR

    # cut streams at surface waves
    rt_cut = rt[0].data[start_sample:end_sample]
    trv_cut = trv_acc[0].data[start_sample:end_sample]

    step = 1
    baz_list = np.linspace(0, int(360 - step), int(360 / step)) # BAz array
    
    # iterate over all BAz values, correlate in time windows
    corr_list = []
    for i in range(len(baz_list)):
        for j in range(len(rt_cut) // rt_TW):
            corr = correlate(a = rt_cut[j*rt_TW:(j+1)*rt_TW],
                             b = trv_cut[j*trv_TW:(j+1)*trv_TW],
                             shift = 0)  

            corr_list.append(corr[0])

    corr_list = np.asarray(corr_list)
    corr_list = corr_list.reshape(
                        len(baz_list),int(round(len(corr_list)/len(baz_list))))
    
    # iterate over correlations, choose those > 0.9
    corrsum_list = []
    for k in range(len(corr_list[:,0])):
        bazsum = []
        for l in range(len(corr_list[0, :])):
            if corr_list[k, l] >= 0.9:
                bazsum.append(corr_list[k, l])
            else:
                bazsum.append(0.0)

        bazsum = np.asarray(bazsum)
        bazsum = sum(bazsum)
        corrsum_list.append(bazsum)

    # best backazimuth estimated, and maximum correlations for that BAz
    corrsum_list = np.asarray(corrsum_list)
    best_ebaz = baz_list[corrsum_list.argmax()] 
    max_ebaz_xcoef = np.max(corr_list[int(best_ebaz)]) 

    # Determine estimated Backazimuth
    if max(corrsum_list) == 0.0:
        EBA = np.nan
    else:
        EBA = best_ebaz

    return corrsum_list, baz_list, max_ebaz_xcoef, EBA


def get_phase_vel(rt, trv_acc, sec, corrcoefs, surf_start, band_check):

    """
    Calculates the phase velocities and the estimated backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type sec: int
    :param sec: Time window length.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type rotate: :class: `~obspy.core.stream.Stream`
    :param rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :type corrsum: list of floats
    :param corrsum: Sum of all correlation coefficients above a certain
        threshold (0.9) within S-waves and surface waves for each backazimuth.
    :type baz_list: numpy.ndarray
    :param baz_list: Vector containing backazimuths (step: 1°).
    :rtype phasv: numpy.ndarray
    :return phasv: Phase velocities of the seismic signal.
    :rtype EBA: float
    :return EBA: Estimated backazimuth.
    """
    # time window in samples
    rt_TW = int(rt[0].stats.sampling_rate * sec)
    trv_TW = int(trv_acc[0].stats.sampling_rate * sec)

    # if dealing with freq bands, start at ind_surf
    if band_check:  
        start = surf_start
    elif not band_check:
        start = 0

    # calculate phase velocity (km/s) for correlations >= 0.75
    phasv_list = []
    for i in range(start, len(corrcoefs)):
        if corrcoefs[i] >= 0.75:
            phasv = (max(trv_acc[0].data[i*trv_TW:(i+1)*trv_TW]) / 
                     max(rt[0].data[i*rt_TW:(i+1)*rt_TW]))
            phasv *= (1/2) * (1E-3)
        else:
            phasv = np.NaN

        phasv_list.append(phasv)

    phasv_list = np.asarray(phasv_list)  

    return phasv_list


def sn_ratio(full_signal, p_arrival, sampling_rate):

    """
    Characterizes the signal-to-noise ratio of the event(s) as the ratio of
    the peak amplitude of the whole wave train and the mean amplitude in a
    noise window before the first theoretical arrival, assuming that the noise
    has the same behavior in all the data. The inputs are the data, the
    theoretical time of the first P-arrival (as seconds after the first sample
    of the fetched data) and the sampling rate.

    :type full_signal: numpy.ndarray
    :param full_signal: Amplitude data of the full signal.
    :type p_arrival: float
    :param p_arrival: Arrival time of the first P-wave.
    :type sam_rate: float
    :param sam_rate: Sampling rate.
    :rtype SNR: float
    :return SNR: Signal-to-noise ratio of the seismogram.
    """
    SR = int(sampling_rate)
    p_arrival = int(round(p_arrival))

    tr_sign = max(full_signal)
    tr_noise = abs(np.mean(full_signal[SR * (p_arrival - 180): 
                                       SR * (p_arrival - 100)]))
    
    SNR = tr_sign/tr_noise

    return SNR


def store_info_json(rt, ac, trv_acc, corrcoefs, dist_baz, arriv_p, EBA, station, 
                    phasv_means, phasv_stds, startev, event, data_sources, 
                    depth, max_ebaz_xcoef, folder_name, tag_name):

    """
    Generates a human readable .json file that stores data for each event,
    like peak values (acceleration, rotation rate,
    zero-lag correlation coefficient), signal-to-noise ratio, backazimuths.

    :type rotate: :class: `~obspy.core.stream.Stream`
    :param rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type baz: tuple
    :param baz: [0] Great circle distance in m, 
                [1] azimuth A->B in degrees,
                [2] backazimuth B->A in degrees.
    :type arriv_p: float
    :param arriv_p: Arrival time of the first P-wave.
    :type EBA: float
    :param EBA: Estimated backazimuth.
    :phasv_means
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :type srcRT: str
    :param station: Data source for RoTations
    :type srcTR: str
    :param srcTR: Data source for TRanslations
    :type folder_name: string
    :param folder_name: Name of the folder containing the event.
    :type tag_name: string
    :param tag_name: Handle of the event.
    """

    orig = event.preferred_origin() or event.origins[0] # Event origin
    catalog = orig.creation_info.author or orig.creation_info.agency_id
    magnitude = event.preferred_magnitude() or event.magnitudes[0] # Mag info.

    PAT = max(trv_acc[0].data)  # Peak transverse acceleration [nm/s]
    PRZ = max(rt[0].data)  # Peak vertical rotation rate [nrad/s]
    PCC = max(corrcoefs)  # Peak correlation coefficient
    MCC = min(corrcoefs) # Minimum correlation coefficient
    TBA = dist_baz[2]  # Theoretical backazimuth [°]
    DS_KM = 0.001 * dist_baz[0] # Epicentral Distance [km]
    DS_DEG =  DS_KM / 111.11 # Epicentral Distance [°]
    SNT = sn_ratio(ac[0].data, arriv_p, ac[0].stats.sampling_rate) 
    SNR = sn_ratio(rt[0].data, arriv_p, rt[0].stats.sampling_rate) 

    # common event dictionary
    dic_event = OrderedDict([
                ('event_id', event.resource_id.id),
                ('event_source', catalog),
                ('event_latitude', orig.latitude),
                ('event_longitude', orig.longitude),
                ('origin_time', str(startev)),
                ('trace_start', str(startev-180)),
                ('trace_end', str(startev+3*3600)),
                ('magnitude', magnitude.mag),
                ('magnitude_type', magnitude.magnitude_type),
                ('depth', depth),
                ('depth_unit', 'km')
                ])

    # individual station dictionary w/ rotational parameters
    dic_station = OrderedDict([
            ('station_information_{}'.format(station), 
                OrderedDict([
                ('station_latitude', rt[0].stats.coordinates.latitude),
                ('station_longitude', rt[0].stats.coordinates.longitude),
                ('rotation_station', 
                    OrderedDict([
                    ('network', rt[0].stats.network),
                    ('station', rt[0].stats.station),
                    ('loc', rt[0].stats.location),
                    ('channel', rt[0].stats.channel),
                    ('data_source', data_sources['BJZ'])
                    ])
                ),
                ('translation_station', 
                    OrderedDict([
                    ('network', trv_acc[0].stats.network),
                    ('station', trv_acc[0].stats.station),
                    ('loc', trv_acc[0].stats.location),
                    ('channel_N', ac.select(component='N')[0].stats.channel),
                    ('channel_E', ac.select(component='E')[0].stats.channel),
                    ('channel_Z', ac.select(component='Z')[0].stats.channel),
                    ('data_source', data_sources['BHN'])
                    ])
                ),
                ('rotational_parameters',
                    OrderedDict([
                    ('epicentral_distance', DS_KM),
                    ('epicentral_distance_unit', 'km'),
                    ('epicentral_distance_in_deg', DS_DEG),
                    ('theoretical_backazimuth', TBA),
                    ('estimated_backazimuth', EBA),
                    ('backazimuth_unit', 'degrees'),
                    ('max_xcoef_for_estimated_backazimuth', max_ebaz_xcoef),
                    ('peak_vertical_rotation_rate', PRZ),
                    ('peak_vertical_rotation_rate_unit', 'nrad/s'),
                    ('peak_transverse_acceleration', PAT),
                    ('peak_transverse_acceleration_unit', 'nm/s^2'),
                    ('peak_correlation_coefficient', PCC),
                    ('minimum_correlation_coefficient', MCC),
                    ('vertical_rotation_rate_SNR', SNR),
                    ('transverse_acceleration_SNR', SNT),
                    ])
                ),
                ('phase_velocities', 
                    OrderedDict([
                    ('band_1', 
                        OrderedDict([
                        ('freqmin', 0.01),
                        ('freqmax', 0.02),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[0]),
                        ('vel_std', phasv_stds[0]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_2', 
                        OrderedDict([
                        ('freqmin', 0.02),
                        ('freqmax', 0.04),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[1]),
                        ('vel_std', phasv_stds[1]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_3', 
                        OrderedDict([
                        ('freqmin', 0.04),
                        ('freqmax', 0.10),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[2]),
                        ('vel_std', phasv_stds[2]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_4', 
                        OrderedDict([
                        ('freqmin', 0.10),
                        ('freqmax', 0.20),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[3]),
                        ('vel_std', phasv_stds[3]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_5', 
                        OrderedDict([
                        ('freqmin', 0.20),
                        ('freqmax', 0.30),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[4]),
                        ('vel_std', phasv_stds[4]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_6', 
                        OrderedDict([
                        ('freqmin', 0.30),
                        ('freqmax', 0.40),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[5]),
                        ('vel_std', phasv_stds[5]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_7', 
                        OrderedDict([
                        ('freqmin', 0.40),
                        ('freqmax', 0.60),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[6]),
                        ('vel_std', phasv_stds[6]),
                        ('vel_unit', 'km/s')
                        ])
                    ),
                    ('band_8', 
                        OrderedDict([
                        ('freqmin', 0.60),
                        ('freqmax', 1.0),
                        ('freq_unit', 'Hz'),
                        ('mean_phase_vel', phasv_means[7]),
                        ('vel_std', phasv_stds[7]),
                        ('vel_unit', 'km/s')
                        ])
                    )
                    ])
                )
                ])
            )
            ])  

    # if: json already created for previous station, overwrite event info
    # else: write a new json file !!! assumes the event information is the same
    filename_json = os.path.join(folder_name,tag_name + '.json')

    if os.path.exists(filename_json): 
        dic_event = json.load(open(filename_json),object_pairs_hook=OrderedDict)

    # combine two dictionaries, save
    dic_event.update(dic_station)

    outfile = open(filename_json, 'wt')
    json.dump(dic_event,outfile,indent=4)
    outfile.close()


def store_info_xml(event,folder_name,tag_name,station):

    """
    Write QuakeML file. Store extra parameters under the namespace rotational
    seismology. Stations are taken care of in nested tags in the extra tag
    of the xml file. Parameters used for filtering events on the JANE database
    framework, taken from .json file 

    :type folder_name: string
    :param folder_name: Name of the folder containing the event.
    :type tag_name: string
    :param tag_name: Handle of the event.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    """
    ns = 'http://www.rotational-seismology.org'
    filename_json = os.path.join(folder_name,tag_name + '.json')
    filename_xml = os.path.join(folder_name,tag_name + '.xml')

    # check if xml already exists
    if os.path.exists(filename_xml):
        event = read_events(filename_xml,format='QUAKEML')[0]

    # check if xml already has extra section
    try:
        event.extra
    except AttributeError:
        event.extra = AttribDict()

    # grab data from json file
    data = json.load(open(filename_json))
    rotational_parameters = ['epicentral_distance',
                             'theoretical_backazimuth',
                             'peak_correlation_coefficient']
    
    # write parameters from json into attribute dictionaries
    params = AttribDict()
    for RP in rotational_parameters:
        RP_value = (data['station_information_{}'.format(station)]
                                                ['rotational_parameters'][RP])
        params[RP] = AttribDict()
        params[RP] = {'namespace':ns,
                      'value': RP_value}

    # set unit attributes
    params.epicentral_distance.attrib = {'unit':"km"}
    params.theoretical_backazimuth.attrib = {'unit':"degree"}

    event.extra['rotational_parameters_{}'.format(station)] = \
                                                            {'namespace': ns,
                                                            'value': params}

    event.write(filename_xml,'QUAKEML',
                        nsmap={"rotational_seismology_database": 
                            r"http://www.rotational-seismology.org"})



def plot_waveform_comp(event, station, link, mode, folder_name, tag_name):

    """
    Compares vertical rotation rate and transversal acceleration through
    direct waveform comparison in different time windows and through cross-
    correlation analysis. It also stores some values obtained through the
    routine, like peak values (signal amplitudes, correlation coefficients)
    and signal-to-noise ratios

    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :type link: string
    :param link: Link to the Iris-xml, where event and moment tensor data are
        fetched.
    :type mode: str
    :param mode: Defines where data fetched from
    :type folder_name: string
    :param folder_name: Name of the folder containing the event.
    :type tag_name: string
    :param tag_name: Handle of the event.
    """

    # gather event information, stream obsjects etc.
    event_lat, event_lon, depth, startev, rt, ac, dist_baz, data_sources = \
        event_info_data(event, station, mode, polarity, instrument)
    
    # distance in km, backazimuth in deg
    DS = dist_baz[0]
    BAz = dist_baz[2]

    station_lat = rt[0].stats.coordinates.latitude
    station_lon = rt[0].stats.coordinates.longitude

    # try to get region information
    try:
        region = event.event_descriptions[0]['text']
        reg = True
    except:
        reg = False
    
    if link != 'blank':
        MomentTensor, Magnitude, Region = get_moment_tensor_magnitude(link)

    # =========================================================================
    #                                   PAGE 1
    #               Create map with event location & information
    #
    # =========================================================================
    print("\nPage 1 > Title Card...", end=" ")

    if is_local(DS) == 'non-local': 
        if DS <= 13000000:
            plt.figure(figsize=(18, 9))
            plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

            # globe plot
            map = Basemap(projection='ortho', 
                            lat_0=(station_lat + event_lat) / 2, 
                            lon_0=(station_lon + event_lon) / 2, 
                            resolution='l')
            map.drawmeridians(np.arange(0, 360, 30))
            map.drawparallels(np.arange(-90, 90, 30))
        elif DS > 13000000:
            plt.figure(figsize=(18, 9))
            plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

            # If the great circle between the station and event is crossing the
            # 180° meridian in the pacific and the stations are far apart the map
            # has to be re-centered, otherwise wrong side of globe.
            if abs(station_lon - event_lon) > 180:
                lon0 = 180 + (station_lon + event_lon) / 2
            else:
                lon0 = (station_lon + event_lon) / 2

            map = Basemap(projection='moll', lon_0=lon0, resolution='l')
            map.drawmeridians(np.arange(0, 360, 30))
            map.drawparallels(np.arange(-90, 90, 30))

    elif is_local(DS) == 'local':
        plt.figure(figsize=(18, 9))
        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

        # conic map plot
        map = Basemap(projection='lcc', 
                        lat_0=(station_lat + event_lat) / 2, 
                        lon_0=(station_lon + event_lon) / 2, 
                        resolution='i', 
                        width=3000000, 
                        height=2000000)
        map.drawparallels(np.arange(0., 90, 5.), labels=[1, 0, 0, 1])
        map.drawmeridians(np.arange(0., 360., 5.), labels=[1, 0, 0, 1])
        map.drawstates(linewidth=0.25)

    elif is_local(DS) == 'close':
        plt.figure(figsize=(26, 13))
        plt.title('Event: %s %sZ\n \n ' %
                (startev.date, startev.time), fontsize=24, fontweight='bold')
        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

        # conic map plot
        map = Basemap(projection='lcc', 
                        lat_0=(station_lat + event_lat) / 2, 
                        lon_0=(station_lon + event_lon) / 2, 
                        resolution='i', 
                        width=600000, 
                        height=400000)
        map.drawparallels(np.arange(0., 90, 2.), labels=[1, 0, 0, 1])
        map.drawmeridians(np.arange(0., 360., 2.), labels=[1, 0, 0, 1])
        map.drawrivers(linewidth=0.25, color='b')
        map.drawstates(linewidth=0.25)

    # basemap boundary settings
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral', lake_color='lightblue')
    map.drawmapboundary(fill_color='lightblue')

    if is_local(DS) == 'local' or is_local(DS) == 'close':
        map.drawlsmask(land_color='coral', ocean_color='lightblue', lakes=True)
        map.drawcountries(linewidth=0.6)
   
    map.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, 
                                                    linewidth=3, color='yellow')

    # add beachballs for the event and station triangle
    if station == 'RLAS':
        x, y = map(event_lon, event_lat)
        statlon, statlat = map(station_lon, station_lat)

        if is_local(DS) == 'non-local':
            map.scatter(statlon, statlat, 200, color="b", marker="v",
                                                     edgecolor="k", zorder=100)
            plt.text(statlon + 200000, statlat, station, va="top",
                                family="monospace", weight="bold", zorder=101,
                                color='k', backgroundcolor='white')
            map.scatter(x, y, 200, color="b", marker="*", edgecolor="k", 
                                                                    zorder=100)

            if mode == 'link':
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s \n %s \n \n' % (startev, Region),
                                                    fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')
                b = Beach(MomentTensor, xy=(0.5, 0.5), facecolor='blue',
                                            width=0.5, linewidth=1, alpha=1.0)
                b.set_zorder(200)
                ax.add_collection(b)

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]'
                          '\n\nDepth: %.2f [km]'
                          % (str(event).split('\n')[0][57:64], 0.001 * DS,
                             0.001 * DS / 111.11, depth), fontsize=18,
                            fontweight='bold')
                ax = plt.gca()
                ax.axis('off')

                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor'
                          ' Catalog (GCMT)' 
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

            else:
                map.scatter(x, y, 200, color="b", marker="*", edgecolor="k",
                            zorder=100)

                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s %sZ\n' % (startev.date, startev.time), 
                                                    fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                if reg==True:
                    plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]\n\nDepth: %.2f [km]'
                              % (region, str(event).split('\n')[0][57:64], 0.001 * DS,
                                 0.001 * DS / 111.11, depth),
                                fontsize=18, fontweight='bold')
                else:
                    plt.title(u'\n\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]\n\nDepth: %.2f [km]'
                              % (str(event).split('\n')[0][57:64], 0.001 * DS,
                                 0.001 * DS / 111.11, depth),
                                fontsize=18, fontweight='bold')
                ax = plt.gca()
                ax.axis('off')

                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor'
                          ' Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

        else:
            map.scatter(statlon, statlat, 200, color="b", marker="v",
                                                    edgecolor="k", zorder=100)
            plt.text(statlon + 35000, statlat, station, fontsize=12, va="top",
                                family="monospace", weight="bold", zorder=101,
                                color='k', backgroundcolor='white')
            map.scatter(x, y, 300, color="b", marker="*", edgecolor="k", 
                                                                    zorder=100)

            if mode == 'link':
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s \n %s \n \n' % (startev, Region),
                          fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')
                b = Beach(MomentTensor, xy=(0.5, 0.5), facecolor='blue',
                                            width=0.5, linewidth=1, alpha=1.0)
                b.set_zorder(200)
                ax.add_collection(b)

                plt.subplot2grid((4, 9), (2, 0), colspan=2)

                plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]'
                          '\n\nDepth: %.2f [km]'
                          % (str(event).split('\n')[0][57:64], 0.001 * DS,
                             0.001 * DS / 111.11, depth), fontsize=18,
                          fontweight='bold')
                ax = plt.gca()
                ax.axis('off')
                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor'
                          ' Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

            else:
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s %sZ\n' % (startev.date, startev.time), 
                                                    fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                if reg==True:
                    plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]\n\nDepth: %.2f [km]'
                              % (region, str(event).split('\n')[0][57:64], 0.001 * DS,
                                 0.001 * DS / 111.11, depth),
                                                 fontsize=18, fontweight='bold')
                else:
                    plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]'
                              '\n\nDepth: %.2f [km]'
                              % (str(event).split('\n')[0][57:64], 0.001 * DS,
                                 0.001 * DS / 111.11, depth), 
                                                fontsize=18,fontweight='bold')
                ax = plt.gca()
                ax.axis('off')
                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor'
                          ' Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

    # stations other than RLAS, not tested yet
    else:
        x, y = map(event_lon, event_lat)
        statlon, statlat = map(station_lon, rt[0].stats.
                               coordinates.latitude)

        if is_local(DS) == 'local':
            map.scatter(x, y, 600, color="b", marker="*", edgecolor="k",
                                                                    zorder=100)
            map.scatter(statlon, statlat, 700, color="b", marker="v",
                                                    edgecolor="k", zorder=100)
            plt.text(statlon + 27000, statlat, station, fontsize=18, va="top",
                                family="monospace", weight="bold", zorder=101,
                                 color='k', backgroundcolor='white')
        elif is_local(DS) == 'non-local':
            map.scatter(x, y, 200, color="b", marker="*", edgecolor="k",
                                                                     zorder=100)
            map.scatter(statlon, statlat, 300, color="b", marker="v",
                                                    edgecolor="k", zorder=100)
            plt.text(statlon + 200000, statlat, station, va="top",
                                family="monospace", weight="bold", zorder=101,
                                color='k', backgroundcolor='white')
        else:
            map.scatter(x, y, 250, color="b", marker="*", edgecolor="k",
                                                                     zorder=100)
            map.scatter(statlon, statlat, 400, color="b", marker="v",
                                                      edgecolor="k", zorder=100)
            plt.text(statlon + 6000, statlat + 1000, station, va="top",
                                  family="monospace", weight="bold", zorder=101,
                             color='k', backgroundcolor='white', fontsize=18)

        plt.subplot2grid((4, 9), (1, 0), colspan=2)
        plt.title(u'Event: %s %sZ\n' % (startev.date, startev.time), 
                                                    fontsize=20, weight='bold')
        ax = plt.gca()
        ax.axis('equal')
        ax.axis('off')

        plt.subplot2grid((4, 9), (3, 0), colspan=2)
        if reg==True:
            plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]\n\nDepth: %.2f [km]'
                      % (region, str(event).split('\n')[0][57:64], 0.001 * DS,
                         0.001 * DS / 111.11, depth),
                      fontsize=18, fontweight='bold')
        else:
            plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km], %.2f [°]'
                      '\n\nDepth: %.2f [km]'
                      % (str(event).split('\n')[0][57:63], 0.001 * DS,
                         0.001 * DS / 111.11, depth), fontsize=18,
                      fontweight='bold')
        ax = plt.gca()
        ax.axis('off')

    # ============================= Save Figure ============================
    plt.savefig(
        os.path.join(folder_name,tag_name + '_{}_page_1.png'.format(station)))
    plt.close()
    print("Done")

    # ======================================================================== 
    #                                   
    #               Preprocessing of rotations and translations
    #
    # =========================================================================
    rt, ac, rt_pcoda, ac_pcoda, sec, sec_p, cutoff, cutoff_pc = resample(
                                                    is_local(DS), rt, ac)

    print("Removing instrument response...")
    rt, ac, rt_pcoda, ac_pcoda = remove_instr_resp(
                                    rt, ac, rt_pcoda, ac_pcoda,station, startev)

    print("Filtering and rotating traces...")
    # filter raw data, rotate some to theoretical backazimuth
    rt_bands, rotate_bands, rotate, rt_pcoda_coarse, pcoda_rotate_coarse, \
    filt_rt_pcoda, filt_ac_pcoda, filt_rotate = filter_and_rotate(ac, rt, BAz, 
                rt_pcoda, ac_pcoda, cutoff, cutoff_pc, station, is_local(DS))

    print("Getting theoretical arrival times...")
    # find event start and theoretical arrival times for seismics phases
    init_sec = startev - ac[0].stats.starttime
    arriv_p, arriv_s = ps_arrival_times(DS, depth, init_sec)
    min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf = \
                    time_windows(DS, arriv_p, arriv_s, init_sec, is_local(DS))

    # ======================================================================== 
    #                                
    #                                   Page 2
    #           Waveform Comparison plot w/ individual phase subplots
    #
    # =========================================================================    
    print("\nPage 2 >  Waveform Comparison...",end=" ")

    # separate transverse acceleration from rotated traces
    trv_acc = rotate.select(component='T')

    # rt.taper(max_percentage=0.05) # this was here but we already taper?
    c1 = .5 * max(abs(trv_acc[0].data)) / max(abs(rt[0].data))  # vel in m/s
    fact1 = 2 * max(rt[0].data)  # Factor to displace the rotation rate
    time = rt[0].stats.delta * np.arange(0, len(rt[0].data))  # time in seconds

    # main figure
    plt.figure(figsize=(18, 9))
    plt.subplot2grid((6, 5), (2, 0), colspan=5, rowspan=2)
    plt.plot(time, rt[0].data, color='r', label=r'Rotation rate')
    plt.plot(time, (0.5 / c1) * trv_acc[0].data + fact1, 
                                color='k', label=r'Transversal acceleration')
    plt.xlabel(r'Time [s]', fontweight='bold', fontsize=13)
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=13)
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(min(rt[0].data), fact1 + max((1. / (2. * c1)) * trv_acc[0].data))
    box_yposition = ((fact1 + max((1. / (2. * c1)) * trv_acc[0].data))
                     - abs(min(rt[0].data)))/2  # box is in middle of figure
    
    if is_local(DS) == 'non-local':  # gap between annotation and vline
        xgap = 50
    else:
        xgap = 15

    bbox_props = dict(boxstyle="square, pad=0.3", fc='white')
    plt.axvline(x=min_pw, linewidth=1)
    plt.annotate('1', xy=(min_pw+xgap, box_yposition), fontsize=14,
                                            fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_sw, linewidth=1)
    plt.annotate('2', xy=(min_sw+xgap, box_yposition), fontsize=14,
                                            fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_lwi, linewidth=1)
    plt.annotate('3', xy=(min_lwi+xgap, box_yposition), fontsize=14,
                                            fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_lwf, linewidth=1)
    plt.annotate('4', xy=(min_lwf+xgap, box_yposition), fontsize=14,
                                            fontweight='bold', bbox=bbox_props)
    plt.title(r'Ring laser and broadband seismometer recordings. Event: %s %sZ'
              % (startev.date, startev.time))
    plt.grid(True)
    plt.legend(loc=7,shadow=True)

    # =============================== P coda ===============================
    plt.subplot2grid((6, 5), (0, 0), colspan=2)
    
    # integer sampling rate for slice indexing
    rt_SR = int(rt[0].stats.sampling_rate)
    rt_pcoda_SR = int(rt_pcoda_coarse[0].stats.sampling_rate)
    
    # transverse acceleration pcoda
    pc_trv_acc = pcoda_rotate_coarse.select(component='T')

    # pwave arrival times in samples for normal trace and pcoda
    min_pw_rt = rt_SR * min_pw
    max_pw_rt = rt_SR * max_pw
    min_pw_pcrt = rt_pcoda_SR * min_pw
    max_pw_pcrt = rt_pcoda_SR * max_pw

    # pcoda phase velocity
    cp = 0.5 * (max(abs(pc_trv_acc[0][min_pw_pcrt:max_pw_pcrt]))/
                max(abs(rt_pcoda_coarse[0].data[min_pw_pcrt:max_pw_pcrt])))

    # find min and max trace amplitudes for y-limits
    min_ta_pcod = min((0.5 / cp) * pc_trv_acc[0][min_pw_rt:max_pw_rt])
    max_ta_pcod = max((0.5 / cp) * pc_trv_acc[0][min_pw_rt:max_pw_rt])
    min_rt_pcod = min(rt_pcoda_coarse[0].data[min_pw_rt:max_pw_rt])
    max_rt_pcod = max(rt_pcoda_coarse[0].data[min_pw_rt:max_pw_rt])

    plt.plot(time, rt_pcoda_coarse[0].data, color='r')
    plt.plot(time, (0.5 / cp) * pc_trv_acc[0], color='k')
    plt.xlim(min_pw, max_pw)
    plt.ylim(min([min_ta_pcod, min_ta_pcod]), max([max_rt_pcod, max_rt_pcod]))
    plt.xlabel(r'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(u'1: P-coda (Highpass, cut-off: %.1f Hz)' % cutoff_pc)
    plt.grid(True)

    # =============================== S wave ===============================
    plt.subplot2grid((6, 5), (0, 3), colspan=2)

    # swave arrival times in samples
    min_sw_rt = rt_SR * min_sw
    max_sw_rt = rt_SR * max_sw

    cs = 0.5 * (max(abs(trv_acc[0].data[min_sw_rt:max_sw_rt])) / 
                max(abs(rt[0].data[min_sw_rt:max_sw_rt])))

    max_rt_s = max(rt[0].data[min_sw_rt:max_sw_rt])
    min_rt_s = min(rt[0].data[min_sw_rt:max_sw_rt])
    min_ta_s = min((0.5 / cs) * trv_acc[0].data[min_sw_rt:max_sw_rt])
    max_ta_s = max((0.5 / cs) * trv_acc[0].data[min_sw_rt:max_sw_rt])

    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cs) * trv_acc[0].data, color='k')
    plt.xlim(min_sw, max_sw)
    plt.ylim(min([min_ta_s, min_rt_s]), max([max_ta_s, max_rt_s]))
    plt.xlabel(r'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(u'2: S-wave (Lowpass, cut-off: %s Hz)' % (cutoff))
    plt.grid(True)

    # ======================== Initial surface waves ==========================
    plt.subplot2grid((6, 5), (5, 0), colspan=2)

    # initial surface wave arrival times in samples
    min_lwi_rt = rt_SR * min_lwi
    max_lwi_rt = rt_SR * max_lwi

    cl1 = 0.5 * (max(abs(trv_acc[0].data[min_lwi_rt:max_lwi_rt])) /
                 max(abs(rt[0].data[min_lwi_rt:max_lwi_rt])))

    min_rt_lwi = min(rt[0].data[min_lwi_rt:max_lwi_rt])
    max_rt_lwi= max(rt[0].data[min_lwi_rt:max_lwi_rt])
    min_ta_lwi = min((0.5 / cl1) * trv_acc[0].data[min_lwi_rt:max_lwi_rt])
    max_ta_lwi = max((0.5 / cl1) * trv_acc[0].data[min_lwi_rt:max_lwi_rt])

    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cl1) * trv_acc[0].data, color='k')
    plt.xlim(min_lwi, max_lwi)
    plt.ylim(min([min_rt_lwi, min_ta_lwi]),max([max_rt_lwi, max_ta_lwi]))
    plt.xlabel(r'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [rad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(r'3: Initial surface waves (Lowpass, cut-off: %s Hz)'
              % (cutoff))
    plt.grid(True)

    # ========================== Later surface waves ==========================
    plt.subplot2grid((6, 5), (5, 3), colspan=2)

    # later surface wave arrival times in samples
    min_lwf_rt = rt_SR * min_lwf
    max_lwf_rt = rt_SR * max_lwf

    cl2 = 0.5 * (max(abs(trv_acc[0].data[min_lwf_rt:max_lwf_rt])) /
                 max(abs(rt[0].data[min_lwf_rt:max_lwf_rt])))
    
    min_rt_lwf = min(rt[0].data[min_lwf_rt:max_lwf_rt])
    max_rt_lwf = max(rt[0].data[min_lwf_rt:max_lwf_rt])
    min_ta_lwf = min((0.5 / cl2) * trv_acc[0].data[min_lwf_rt:max_lwf_rt])
    max_ta_lwf = max((0.5 / cl2) * trv_acc[0].data[min_lwf_rt:max_lwf_rt])

    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cl2) * trv_acc[0].data, color='k')
    plt.xlim(min_lwf, max_lwf)
    plt.ylim(min([min_rt_lwf, min_ta_lwf]),max([max_rt_lwf, max_ta_lwf]))
    plt.xlabel(r'Time [s]', fontsize=11, fontweight='bold')
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [rad/s] - a$_\mathbf{T}$/2c'
                                        '[1/s]', fontsize=11, fontweight='bold')
    plt.title(r'4: Later surface waves (Lowpass, cut-off: %s Hz)' % (cutoff))
    plt.grid(True)

    # ============================= Save Figure =============================
    plt.savefig(
        os.path.join(folder_name,tag_name + '_{}_page_2.png'.format(station)))
    plt.close()
    print("Done")

    # ======================================================================== 
    #                                
    #                Cross Correlations and Phase Velocities
    #
    # ========================================================================= 
    print("Finding zero-lag correlation coefficients...")

    # correlate vertical rotation rate and transverse acceleration
    corrcoefs, thres = get_corrcoefs(rt, trv_acc, sec, station)

    # calculate correlations for different frequency bands given by seconds_list
    corrcoefs_bands, thresholds = [], []
    seconds_list = [200, 100, 50, 20, 12, 10, 8, 6]

    for i in range(len(rt_bands)):
        tr_acc_band = rotate_bands[i].select(component='T')
        corrcoefs_tmp, thresh_tmp = get_corrcoefs(streamA = rt_bands[i],
                                                  streamB = tr_acc_band,
                                                  sec = seconds_list[i], 
                                                  station = station)
        corrcoefs_bands.append(corrcoefs_tmp)
        thresholds.append(thresh_tmp)

    # zero-lag correlation coefficients for range of backazimuths
    print("Analyzing correlation by BAz bins...")
    corrbaz, maxcorr, backas, max_coefs_10deg = baz_analysis(
                                                    rt, ac, sec, station)

    # estimate backazimuth and correlations for given BAz
    print("Estimating best backazimuth value...")
    corrsum, baz_list, max_ebaz_xcoef, EBA = estimate_baz(
                                        rt, trv_acc, min_sw, max_lwf, station)

    print("Calculating phase velocities...")

    # calculate phase velocities starting at surface wave arrivals 
    surf_start = min_lwi // sec
    phasv = get_phase_vel(
                    rt, trv_acc, sec, corrcoefs, surf_start, band_check=False)
    
    # calculate phase velocities for different frequency bands, 1-8
    phasv_bands,phasv_means,phasv_stds = [],[],[]
    for i in range(len(rt_bands)):
        trv_tmp = rotate_bands[i].select(component='T')
        phasv_tmp = get_phase_vel(rt_bands[i], trv_tmp, seconds_list[i],
                                corrcoefs_bands[i], surf_start, band_check=True)
        
        # filter out NaNs and append to list
        phasv_bands.append(phasv_tmp[~np.isnan(phasv_tmp)])
    
    # phase velocity mean values and std. for json file
    for PVB in phasv_bands:
        if len(PVB) != 0:
            phasv_means.append(np.mean(PVB))
            phasv_stds.append(np.std(PVB))
        else:
            phasv_means.append(np.NaN)
            phasv_stds.append(np.NaN)

    import pdb;pdb.set_trace()
    # ======================================================================== 
    #                                
    #                                Page 3
    #              Cross Correlation, phase velocity determination,
    #                    Estimation of backazimuth figures
    #
    # ========================================================================= 
    print("\nPage 3 > Cross-correlation, Phase velocity...",end=" ")

    X, Y = np.meshgrid(np.arange(0, sec * len(corrcoefs), sec), backas)

    # subfigure 1
    plt.figure(figsize=(18, 9))
    plt.subplot2grid((4, 26), (0, 0), colspan=25)
    plt.plot(time, rt[0].data, color='r', label=r'Rotation rate')
    plt.plot(time, (1. / (2. * c1)) * trv_acc[0].data + fact1, 
                                color='k', label=r'Transversal acceleration')
    plt.ylabel(
        r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c '
        '[1/s]', fontsize=10, fontweight='bold')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(min(rt[0].data), fact1 + max((1. / (2. * c1)) * trv_acc[0].data))
    plt.title(r'Cross-correlation for $\dot\Omega_z$ and a$_T$ in %s seconds '
                          'time windows (lowpass, cutoff: %s Hz). Event: %s %sZ'
                          % (sec, cutoff, startev.date, startev.time))
    plt.grid(True)
    plt.legend(loc=7,shadow=True)

    # subfigure 2
    plt.subplot2grid((4, 26), (1, 0), colspan=25)
    plt.scatter(np.arange(0, sec * len(phasv), sec), phasv,
                c=corrcoefs, vmin=0.75, vmax=1, s=35,
                cmap=plt.cm.autumn_r)
    plt.ylabel(r'Phase velocity [km/s]', fontsize=10, fontweight='bold')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(0, 16)
    plt.grid(True)

    fig = plt.subplot2grid((4, 26), (1, 25))
    cmap = mpl.cm.autumn_r
    norm = mpl.colors.Normalize(vmin=0.75, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(
                            fig, cmap=cmap, norm=norm, orientation='vertical',
                            ticks=[0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')

    # subfigure 3
    plt.subplot2grid((4, 26), (2, 0), colspan=25)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), max_coefs_10deg, 'ro-', 
                                    label='Max. CC for est. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), corrcoefs, 'ko-', 
                                        label='CC for theo. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), thres, '--r', lw=2)
    plt.ylabel(r'X-corr. coeff.', fontsize=10, fontweight='bold')
    plt.text(time[len(time) - 1] + 50, 0.71, r'0.75', color='red')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))

    min_corr = min(min(max_coefs_10deg), min(corrcoefs))
    plt.ylim(min_corr, 1)
    plt.legend(loc=4, shadow=True)
    plt.grid(True)

    # subfigure 4
    plt.subplot2grid((4, 26), (3, 0), colspan=25)
    teobaz = BAz * np.ones(len(corrcoefs) + 1)
    plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r, vmin=-1, vmax=1)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), teobaz, '--r', lw=2)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), maxcorr, '.k')
    plt.text(1000, BAz, str(BAz)[0:5] + r'°',
             bbox={'facecolor': 'black', 'alpha': 0.8}, color='r')

    if EBA == baz_list[np.asarray(corrsum).argmax()]:
        obsbaz = EBA * np.ones(len(corrcoefs) + 1)
        plt.text(400, EBA, str(EBA)[0:5] + r'°',
                 bbox={'facecolor': 'black', 'alpha': 0.8}, color='y')
        plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec),
                 obsbaz, '--y', lw=2)

    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.xlabel(r'Time [s]', fontweight='bold')
    plt.ylabel(r'BAz [°]', fontsize=10, fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks([0,60,120,180,240,300,360])
    plt.grid(True)
    fig = plt.subplot2grid((4, 26), (3, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(fig, cmap=plt.cm.RdYlGn_r, norm=norm, 
                                                        orientation='vertical')
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks([-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0])

    # ============================= Save Figure =============================
    
    plt.savefig(
        os.path.join(folder_name,tag_name + '_{}_page_3.png'.format(station)))
    plt.close()
    print("Done")

    # ======================================================================== 
    #                                
    #                                P-Coda analysis
    #
    # ========================================================================= 
    
    print("Analyzing rotations in the P-coda")
    
    # Zero-lag correlation coefficients
    print("ZLCC...",end=" ")

    # integer sampling rates
    rt_pc_SR = int(filt_rt_pcoda[0].stats.sampling_rate)
    trv_pc_SR = int(filt_rotate[0].stats.sampling_rate)
    
    corrcoefs_p = []
    lwi_average = int(round((min_lwi+max_lwi)/2))

    # separate transverse and vertical components
    trv_acc_pcoda = filt_rotate.select(component='T')
    vrt_acc_pcoda = filt_rotate.select(component='Z')

    # cut pcoda at correct time window and taper cuts
    rt_pcoda_cut = filt_rt_pcoda.copy()
    trv_pcoda_cut = trv_acc_pcoda.copy()
    rt_pcoda_cut[0].data = filt_rt_pcoda[0].data[0:lwi_average * rt_pc_SR]
    trv_pcoda_cut[0].data = trv_pcoda_cut[0].data[0:lwi_average * trv_pc_SR]
    for traces in [rt_pcoda_cut,trv_pcoda_cut]:
        traces.taper(max_percentage=0.05)

    # find correlations
    corrcoefs_p, thres_p = get_corrcoefs(rt_pcoda_cut, trv_pcoda_cut, sec_p,
                                                                        station)

    print("BAz...")

    # surface wave start sample
    max_lwi_ac = trv_pc_SR * max_lwi

    # analyze backazimuth
    corrbaz_p, maxcorr_p, backas, max_coefs_10deg_p = baz_analysis(
                                    rt_pcoda_cut, trv_pcoda_cut, sec_p, station)

    Xp, Yp = np.meshgrid(np.arange(0, sec_p * len(corrcoefs_p), sec_p), backas)

    # set up arrays for plotting
    time_p = rt_pcoda_cut[0].stats.delta * np.arange(0, len(rt_pcoda[0].data))
    fact1_p = 2 * max(rt_pcoda_cut[0].data[0:max_lwi_ac])

    # calculate phase velocity for scaling, from full waveform sections 
    # filtered with the pcoda highpass cutoff
    c1_p = .5 * (max(abs(pc_trv_acc[0].data[0:max_lwi_ac])) /
                 max(abs(rt_pcoda[0].data[0:max_lwi_ac])))

    maxcorr_p_list = []
    for m in range(0, len(maxcorr_p)):
        if np.max(corrbaz_p[:,m]) >= 0.5:
            maxcorr_p_list.append(maxcorr_p[m])
        else:
            maxcorr_p_list.append(0)

    # ======================================================================== 
    #                                
    #                                   Page 4
    #               Cross correlations for P-Coda time window
    #
    # ========================================================================= 

    # P-coda limits
    xlim1 = rt[0].stats.delta * len(rt[0].data) 

    print("\nPage 4 > P-Coda Waveform Comparison...",end=" ")
    plt.figure(figsize=(18, 9))
    plt.subplot2grid((5, 26), (0, 0), colspan=25)
    plt.plot(time_p, vrt_acc_pcoda[0].data, color='g')
    plt.ylabel(r'a$_\mathbf{Z}$ [nm/s$^2$]', fontweight='bold', fontsize=11)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(min(vrt_acc_pcoda[0].data[0:max_lwi_ac]),
             max(vrt_acc_pcoda[0].data[0:max_lwi_ac]))
    plt.title(r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ and a$_\mathbf{T}$'
              'correlation in the P-coda in a %d seconds time window'
              ' (highpass, cutoff: 1 Hz). Event: %s %sZ' % 
                                            (sec_p, startev.date, startev.time))
    plt.axvline(x=min_pw, linewidth=1)
    plt.axvline(x=min_sw, linewidth=1)
    plt.subplot2grid((5, 26), (1, 0), colspan=25, rowspan=2)

    plt.plot(time_p, rt_pcoda[0].data, color='r', label=r'Rotation rate')
    plt.plot(time_p, (0.5 / c1_p) * trv_acc_pcoda[0].data + fact1_p, color='k',
                                             label=r'Transversal acceleration')
    plt.ylabel(r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] -'
                   'a$_\mathbf{T}$/2c [1/s]', fontweight='bold', fontsize=11)
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(min(rt_pcoda[0].data[0:max_lwi_ac]), 
        fact1_p + max((1. / (2. * c1_p)) * trv_acc_pcoda[0].data[0:max_lwi_ac]))
    xlim2 = (min_lwi + max_lwi) // 2
    box_yposition2 = (fact1_p + max((1. / (2. * c1_p)) * 
                        trv_acc_pcoda[0].data[0:max_lwi_ac]) -
                        np.abs(min(rt_pcoda[0].data[0: max_lwi_ac])))/2.
    plt.axvline(x=min_pw, linewidth=1)
    plt.annotate('P-arrival', 
                            xy=(min_pw+xgap*float(xlim2/xlim1),box_yposition2), 
                            fontsize=14, fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_sw, linewidth=1)
    plt.annotate('S-arrival', 
                            xy=(min_sw+xgap*float(xlim2/xlim1),box_yposition2), 
                            fontsize=14,fontweight='bold', bbox=bbox_props)
    plt.grid(True)
    plt.legend(loc=6, shadow=True)


    plt.subplot2grid((5, 26), (3, 0), colspan=25)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p), corrcoefs_p, '.k')
    plt.ylabel(r'X-corr. coeff.', fontweight='bold')
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(0, 1)
    plt.grid(True)
    plt.subplot2grid((5, 26), (4, 0), colspan=25)
    plt.pcolor(Xp, Yp, corrbaz_p, cmap=plt.cm.RdYlGn_r, vmin=-1, vmax=1)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p), 
                                                        maxcorr_p_list, '.k')
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.xlabel(r'Time [s]', fontweight='bold')
    plt.ylabel(r'BAz [°]', fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks([0,60,120,180,240,300,360])
    plt.grid(True)

    fig = plt.subplot2grid((5, 26), (4, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(fig, cmap=plt.cm.RdYlGn_r, norm=norm, 
                                                        orientation='vertical')
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks([-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0])
   
    # ============================= Save Figure ============================== 

    plt.savefig(
        os.path.join(folder_name,tag_name + '_{}_page_4.png'.format(station)))
    plt.close()
    print("Done")

    # ======================================================================== 
    #                                
    #                        Store in Json and XML files
    #
    # ========================================================================= 
    print("\n>> Storing event information in JSON and XML files...",end=" ")
    
    store_info_json(rt, ac, trv_acc, corrcoefs, dist_baz, arriv_p, EBA, station, 
                    phasv_means, phasv_stds, startev, event, data_sources, 
                    depth, max_ebaz_xcoef, folder_name, tag_name)

    store_info_xml(event,folder_name,tag_name,station)

    print("Done\n")


def generate_tags(event):

    """
    Generates all naming schema tags for an event and prints dialog as it does
    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :rtype tag_name: str
    :return tag_name: event tag i.e. 'GCMT_2017-09-23T125302_6.05_OAXACA_MEXICO'
    :rtype folder_name: str
    :return folder_name: folder name tag
    :rtype check_folder_exists: list
    :return check_folder_exists: glob list with identical file names if event 
                                was already processed
    """

    event_information = str(event).split('\n')[0][7:]
    flinn_engdahl = event.event_descriptions[0]['text'].upper()
    print('{}\n{}\n{}\n{}'.format(
                            bars,flinn_engdahl,event_information,bars))

    # create tags for standard filenaming
    # magnitude, always keep at 3 characters long, fill w/ 0's if not
    mag_tag = '{:0^4}'.format(str(event.magnitudes[0]['mag']))

    # Flinn Engdahl region, i.e. SOUTHEAST_OF_HONSHU_JAPAN
    substitutions = [(', ', '_'),(' ','_')]
    for search, replace in substitutions:
        flinn_engdahl = flinn_engdahl.replace(search,replace)

    # remove '.' from end of region name if necessary (i.e. _P.N.G. > _.P.N.G)
    if flinn_engdahl[-1] == '.':
        flinn_engdahl = flinn_engdahl[:-1]

    # ISO861 Time Format, i.e. '2017-09-23T125302Z'
    orig = event.preferred_origin() or event.origins[0]
    time_tag = orig['time'].isoformat()[:19].replace(':','')+'Z'

    # i.e. 'GCMT_2017-09-23T125302_6.05_OAXACA_MEXICO'
    tag_name = '_'.join((catalog,time_tag,mag_tag,flinn_engdahl))

    # i.e. './OUTPUT/GCMT_2017-09-23T125302_6.05_OAXACA_MEXICO/
    folder_name = os.path.join(output_path,tag_name)

    # short tags used to check if an event with the same time tag has 
    # been processed because different catalogs publish diff. magnitudes
    tag_name_short = '_'.join((catalog,time_tag))
    folder_name_short = os.path.join(output_path,tag_name_short)
    check_folder_exists = glob.glob(folder_name_short + '*')

    return tag_name, folder_name, check_folder_exists


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Comparison of transvere\
        acceleration and vertical rotation rate through direct waveform\
        comparison in different time windows, and cross-correlation analysis.')
    parser.add_argument('--station', help='Choice of station: RLAS, ROMY\
        (default is RLAS)', type=str, default='RLAS')
    parser.add_argument('--mode', help='Choose catalog to download events: \
        GCMT catalog for the most up to date catalog. Link for plotting single \
        event beachballs (requires --link to event). ISC QuakeML file for \
        catalog of local/regional events. IRIS for most stable solutions, \
        though recent events might not be present \
        (default: gcmt, else: iscquakeml, iris, link)', type=str,default='GCMT')
    parser.add_argument('--link',help='URL of moment tensor link, to plot \
        beachballs for a single event --mode',type=str,default='blank')
    parser.add_argument('--polarity', help='Flip polarity of rotation data to \
        fix data errors, to be used in specific time windows of catalog rerun \
        (default: normal, otherwise: reverse)',type=str, default='normal')
    parser.add_argument('--instrument', help='Choose instrument if using RLAS,\
        STS2 GPS went down for specific time, so nearby WETR can be used, \
        though data is lower quality than STS2\
        (default: sts2, otherwise: lennartz)',type=str, default='sts2')
    parser.add_argument('--min_magnitude', help='Minimum magnitude for \
        events (default is 3).', type=float or int, default=4.0)
    parser.add_argument('--max_magnitude', help='Maximum magnitude for \
        events (default is 10).', type=float or int, default=10.0)
    parser.add_argument('--min_depth', help='Minimum depth for events in km \
        (default is 0 km). Positive down for IRIS.', type=float or int, default=0.0)
    parser.add_argument('--max_depth', help='Maximum depth for events in km \
        (default is 1000 km for IRIS).', type=float or int, default=1000.0)
    parser.add_argument('--min_latitude', help='Minimum latitude for events.\
        Format +/- 90 decimal degrees (default is -90°).', type=float or int,
                        default=-90.0)
    parser.add_argument('--max_latitude', help='Maximum latitude for events \
        (default is 90°).', type=float, default=90.0)
    parser.add_argument('--min_longitude', help='Minimum longitude for \
        events. Format +/- 180 decimal degrees (default is -180°).',
                        type=float or int, default=-180.0)
    parser.add_argument('--max_longitude', help='Maximum longitude for \
        events (default is 180°).', type=float or int, default=180.0)
    parser.add_argument('--min_datetime', help='Earliest date and time for \
        the search. Format is UTC: yyyy-mm-dd-[T hh:mm:ss]. \
        Example: 2010-02-27T05:00', type=str, default=str(
                        datetime.datetime.now()-datetime.timedelta(hours=168)))
    parser.add_argument('--max_datetime', help='Latest date and time for \
        the search (default is today).',type=str, default=str(
                                                    datetime.datetime.now()))

    args = parser.parse_args()
    station = args.station
    mode = args.mode.upper()
    polarity = args.polarity.lower()
    instrument = args.instrument.upper()
    link = args.link

    # --link should lead to the quakeml webpage of a certain event on IRIS.
    # The fetched xml provides a moment tensor data for the beachballs.
    # i.e. 'http://www.iris.edu/spudservice/momenttensor/736631/quakeml'

    # [default]: get event catalog from GCMT NEW QUICK,
    # or link mode if moment tensor link is given
    if mode == 'GCMT' or (mode =='link' and link != 'blank'):
        catalog = 'GCMT'
        event_source = 'GCMT'
        if UTCDateTime(args.min_datetime) < UTCDateTime(2014,1,1):
            print("\nDownloading events from NDK catalog")
            cat_all = read_events(
                            './populate_database/NDK_events_before2014.ndk')
        else:
            print("\nDownloading events from GCMT NEW QUICK")
            # a link to quick solutions for past year
            cat_all = read_events('http://www.ldeo.columbia.edu/~gcmt/projects/'
                                            'CMT/catalog/NEW_QUICK/qcmt.ndk')

        cat = cat_all.filter('time > '+str(args.min_datetime), 
                            'time < '+str(args.max_datetime),
                            'magnitude >= '+str(args.min_magnitude), 
                            'magnitude <= '+str(args.max_magnitude),
                            # 'depth <= '+str(args.min_depth), 
                            # 'depth >= '+str(args.max_depth),
                            'longitude >= '+str(args.min_longitude), 
                            'longitude <= '+str(args.max_longitude),
                            'latitude >= '+str(args.min_latitude), 
                            'latitude <= '+str(args.max_latitude))

    # get event catalog from IRIS through FDSN webservice
    elif mode == 'IRIS':
        print("\nDownloading events from IRIS")
        catalog = 'GCMT'
        event_source = 'IRIS'
        c = fdsnClient(event_source)
        cat = c.get_events(minmagnitude=args.min_magnitude,
                           maxmagnitude=args.max_magnitude,
                           mindepth=args.min_depth, 
                           # magnitudetype='Mw',
                           maxdepth=args.max_depth,
                           minlatitude=args.min_latitude,
                           maxlatitude=args.max_latitude,
                           minlongitude=args.min_longitude,
                           maxlongitude=args.max_longitude,
                           starttime=args.min_datetime,
                           endtime=args.max_datetime,
                           catalog=catalog)

    elif mode == 'ISCQUAKEML':
        print("\nDownloading events from ISC QuakeML catalog")
        quakeml = './populate_database/extra_events.xml'
        cat = read_events(quakeml, format='QUAKEML')
        event_source = 'ISC'
        catalog = 'ISC'

    else:
        sys.exit('Invalid mode')

    # set file output path
    output_path = './OUTPUT/'
    if not os.path.exists(output_path): 
        os.makedirs(output_path)

    print("%i event(s) downloaded, beginning processing..." % len(cat))
    success_counter = fail_counter = already_processed = 0
    bars = '='*79
    error_list = []
    for event in cat:
        try:
            tag_name, folder_name, check_folder_exists = generate_tags(event)

            # check if current event folder exists
            if check_folder_exists:
                # check if event source is the same, assumes 0 or 1 files found
                if (os.path.basename(check_folder_exists[0]) != 
                                                os.path.basename(folder_name)):
                    print('This event was processed with another mode\n')
                    already_processed += 1
                    continue

                # if new station, run waveform compare again
                try:
                    filename_json = os.path.join(folder_name,tag_name + '.json')
                    data = json.load(open(filename_json))
                    if data['station_information_{}'.format(station)]:
                        print("This event was already processed\n")
                        already_processed += 1
                    else:
                        try:
                            plot_waveform_comp(event, station, link, mode,
                                                        folder_name, tag_name)
                            success_counter += 1

                        # if any error, remove folder, continue
                        except Exception as e:
                            fail_counter += 1
                            print(e)
                            print("Removing incomplete folder...\n")
                            error_list.append(tag_name)
                            shutil.rmtree(folder_name)

                        # if keyboard interrupt, remove folder, quit
                        except KeyboardInterrupt:
                            print("Removing incomplete folder...\n")
                            shutil.rmtree(folder_name)
                            sys.exit()

                # if json not found, folder is incomplete, continue
                except FileNotFoundError:
                    fail_counter += 1 
                    error_list.append(tag_name)
                    print("Incomplete folder found\n")

            
            # event encountered for the first time, create folder, xml, process
            elif not check_folder_exists:  
                os.makedirs(str(folder_name))
                
                # run processing function
                try:
                    plot_waveform_comp(event, station, link, mode,
                                                        folder_name, tag_name)
                    success_counter += 1
                
                # if any error, remove folder, continue
                except Exception as e:
                    fail_counter += 1
                    print(e)
                    print("Removing incomplete folder...\n")
                    error_list.append(tag_name)
                    shutil.rmtree(folder_name)
               
                # if keyboard interrupt, remove folder, quit
                except KeyboardInterrupt:
                    fail_counter += 1
                    print("Removing incomplete folder...\n")
                    shutil.rmtree(folder_name)
                    sys.exit()

        # if error creating tags, continue
        except Exception as e:
            fail_counter += 1
            print("Error in folder/tag name creation; ",e)

    # print end message
    print('{}\n'.format('_'*79))
    print("Catalog complete, no more events to show")
    print("From a total of %i event(s):\n %i was/were successfully processed"
          "\n %i could not be processed \n %i already processed\n" % (
              len(cat), success_counter, fail_counter, already_processed))
    
    # write error log to see events failed, named by search timeframe
    if len(error_list) > 0:
        if not os.path.exists('./errorlogs'):
            os.makedirs('./errorlogs')

        errorlog_name = './errorlogs/wavComp_{}_{}.txt'.format(
                                args.min_datetime[:10],args.max_datetime[:10])
        print("Writing error log: {}".format(errorlog_name))
        with open(errorlog_name,'w') as f:
            f.write("Error Log Created {}\n".format(datetime.datetime.now()))
            for i in error_list:
                f.write('{}\n'.format(i))



# DEBUGGER
# import pdb; pdb.set_trace()
