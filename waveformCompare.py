#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ROTATIONAL SEISMOLOGY ROUTINES. 
Following the theory and assuming atransversely polarized plane wave, 
this script compares the transverse acceleration and vertical rotation 
rate of an event through:

1) direct waveform comparison for different seismic phases (P-coda, S-waves and
surface waves) between transverse acceleration from a broadband 3-component
translation sensor, and vertical rotation rate measured by rotation instruments

2) zero-lag cross-correlation analysis of the waveforms using the theoretical
backazimuth (BAz) of the event. 

3) calculation of local horizontal phase velocity for small time windows,
which show correlation coefficients larger than 0.75. 

4) estimate of BAz by taking correlations at every value of 0 < BAz < 360 

The routine generates as OUTPUT:

+ a QuakeML (.xml) file that stores data for each event, and contains an extra 
tag with rotational parameters, to be used in the JANE format database

+ a human readable (.json) ordered dictionary text file that contains both
event and station information, as well as output data results from processing

+ 4 images showing waveform compariosn, correlations, phase velocities and 
analysis of rotations in he P-Coda time window

Extra Information:

+ quakeML input files can be read directly. Useful for pulling events from other 
catalogs (ISC, ...), Accomplished by setting: --mode iscquakeml

+ default --mode fetches data from the GCMT New Quick catalog or a local .NDK

+ events can be fetched from IRIS FDSN webservice, which is usually faster
and more abundant especially for older events. Set: --mode iris, magnitudes
are MW however, and MWC from GCMT default mode.

+ non-local events bandstoppped for secondary microseism frequencies (5-12s) 

+ only 0.5 hour recordings shown for local and close events

Running The Script:

# fetch data from available sources (archive if working internally or FDSN 
# webservice if external) for the past week, process and produce output
>>> python waveformCompare.py 

# call arguments to choose event parameters such as magnitude and time
>>> python waveformCompare.py --min_magnitude 6 --max_magnitude 7

# i.e. to get the M9.1 2011 Tohoku-Oki earthquake
>>>python waveformCompare.py --min_magnitude 9 --min_datetime 2011-01-01T00:00 \
--mode iris

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
from pprint import pprint
mpl.use('Agg')
from obspy.core import read
import matplotlib.pylab as plt
from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy import read_events, Catalog
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
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
    Downloads channel data from for the desired event day and origin time.
    Parses the instrument ID to find the correct channel to grab from.
    Cascading search: will first check to see if data is available on file 
    (LMU/FFB archives - assumes these paths do not exist on personal computers);
    if paths are not found, resorts to querying FDSN webservice from the list of
    sources given. Returns the data stream, and the source from which it
    was fetched.

    :type origin_time: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param origin_time: event origin time.
    :type instrument_id: str
    :param net: FDSN SEED name for channel (i.e. `BW.RLAS..BJZ`)
    :type soure: list of str's
    :param source: list of URL's of FDSN webservices.
    :rtype st: :class: `~obspy.core.stream.Stream`
    :return st: data fetched as a stream object.    
    :rtype data_source: str
    :return data_source: Source where data was fetched successfully 
    """

    # check paths to see if running on FFB, LMU or neither
    st = None
    dataDir_get = '/bay200/mseed_online/archive/' #FFB
    if not os.path.exists(dataDir_get):
        dataDir_get = '/import/netapp-m-02-bay200/mseed_online/archive/'#LMU            
    if not os.path.exists(dataDir_get):
        dataDir_get = None
    
    net, sta, loc, cha = instrument_id.split('.')
    
    # if data path exists, create file path names for day and day+1
    if dataDir_get:
        print("Fetching {} data from file".format(net))
        fileName = '.'.join((instrument_id,'D',origin_time.strftime('%Y.%j')))
        filePath = os.path.join(dataDir_get, origin_time.strftime('%Y'),
                                net, sta, cha + '.D', fileName)
        
        origin_time2 = origin_time + 86400
        fileName2 = '.'.join((instrument_id,'D',origin_time2.strftime('%Y.%j')))
        filePath2 = os.path.join(dataDir_get, origin_time2.strftime('%Y'),
                                net, sta, chan + '.D', fileName2)

        # if full path exists, read in data, check if time extends to day+1
        if os.path.isfile(filePath):
            data_source = 'Archive'
            if origin_time.hour > 21:
                st = Stream()
                st.extend(read(pathname_or_url = filePath, 
                               starttime = origin_time - 180,
                               endtime = origin_time + 3 * 3600))
                st.extend(read(pathname_or_url = filePath2, 
                               starttime = UTCDateTime(origin_time2.year,
                                                       origin_time2.month,  
                                                       origin_time2.day, 0, 0),
                               endtime = origin_time + 3 * 3600))
                st.merge(method=-1)
            else:
                st = read(pathname_or_url = filePath, 
                          starttime = origin_time - 180,
                          endtime = origin_time + 3 * 3600)    
        else:
            print("\tFile not found: \n\t {} \n".format(filePath))    
    
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
        raise RotationalProcessingException("Data not available for this event")

    # trim full waveform around event
    st.trim(starttime=origin_time-180, endtime=origin_time+3*3600)
    print("\tDownload of {!s} {!s} data successful".format(
              st[0].stats.station, st[0].stats.channel))

    return st, data_source


def event_info_data(event, station, polarity, instrument):

    """
    Extracts event information and generates necessary processing variables
    Calls the download_data() function in order to grab waveforms.
    Assigns correct station information to all stream objects.
    Calculates great circle distance and backazimuth using lat/lon pairs.

    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :type station: str
    :param station: Station to fetch data from.
    :type polarity: str
    :param polarity: ['normal'] or 'reverse' for flipped rotation signals
    :type instrument: str
    :param instrument: 'WET' or 'WETR' choice for comparison to 'RLAS'
    :rtype event_lat: float
    :return event_lat: Latitude of event in degrees.
    :rtype event_lon: float
    :return event_lon: Longitude of event in degrees.
    :rtype depth: float
    :return depth: Hypocenter depth in km
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :return startev: Event origin time.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Three component broadband translation in N/E/Z components.
    :rtype dist_baz: tuple
    :return dist_baz: [0] great circle distance in m, 
                 [1] theoretical azimuth,
                 [2] theoretical backazimuth.
    :rtype data_sources: dictionary
    :return data_sources: collection of data source for each channel
    """

    # find event start
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
        
    # set attributes necessary for all stations
    event_lat = origin.latitude
    event_lon = origin.longitude
    depth = origin.depth * 0.001  # Depth in km
    dist_baz = gps2dist_azimuth(lat1 = event_lat, lon1 =  event_lon, 
                                lat2 = station_lat, lon2 = station_lon)

    for ca in [ac[0], ac[1], ac[2], rt[0]]:
        ca.stats.coordinates = AttribDict()
        ca.stats.coordinates['longitude'] = station_lon
        ca.stats.coordinates['latitude'] = station_lat
        ca.stats['back_azimuth'] = dist_baz[2]
        ca.stats['starttime'] = startev - 180
        ca.stats['sampling_rate'] = 20.
    
    return event_lat, event_lon, depth, startev, rt, ac, dist_baz, data_sources


def is_local(ds_in_km):

    """
    Check whether the event is close, local or far. For epicentral distance x,
    close:  0° < x < 3°     or        0 km < x < 333.33 km
    local: 3° <= x < 10°    or     333.33 <= x < 1111.1 km
    far:  10° <= x          or  1111.1 km <= x
    
    :type ds_in_km: float
    :param ds_in_km: Event-station distance in km
    :rtype: str
    :return: Self-explaining string for event distance.
    """
    # approximation of distance in degrees
    distance_in_deg = ds_in_km / 111.11

    if distance_in_deg < 10.0:
        if distance_in_deg < 3.0:
            is_local = 'CLOSE'
        else:
            is_local = 'LOCAL'
    else:
        is_local = 'FAR'

    return is_local


def get_moment_tensor(event):

    """
    Extract moment tensor from event class if mode GCMT, to plot beachballs

    :type event: :class: `~obspy.core.event.event.Event`
    :param event: Event information container
    :rtype moment_tensor: list of floats
    :return moment_tensor: List of the six independent MT components
    """
    foc_mech = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    tensor = foc_mech.moment_tensor.tensor
    moment_tensor = []
    for component in ['m_rr','m_tt','m_pp','m_rt','m_rp','m_tp']:
        moment_tensor.append(tensor[component])

    return moment_tensor


def resample(is_local, rt, ac):

    """
    Resample signal dependent on locality of event.

    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband translation
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Resampled rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Resampled three component broadband station signal.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: (Decimated) copy of rt for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: (Decimated) copy of ac for p-coda calculation.
    :rtype sec/sec_p: int
    :return sec/sec_p: Time window length.
    :rtype cutoff/cutoff_pc: float
    :return cutoff: Cut-off frequency for lowpass filter.
    """

    cutoff_pc = 0.5 # cutoff for pcoda lowpass
    if is_local == 'FAR':
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt_pcoda.decimate(factor=2)
        ac_pcoda.decimate(factor=2)
        rt.decimate(factor=4)
        ac.decimate(factor=4)
        sec = 120 # length of time window in seconds
        sec_p = 5 # pcoda time window in seconds
        cutoff = 1.0 # cut-off freq for full-trace lowpass
    elif is_local == 'LOCAL':
        for trr in (rt + ac):
            trr.data = trr.data[0: int(1800 * rt[0].stats.sampling_rate)]
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 5 
        sec_p = 2
        cutoff = 2.0  
    elif is_local == 'CLOSE':
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
    Remove instrument response from the original signal
    General highpass and lowpass to remove signal outside instrument capability
    * Note for STS2 Poles and Zeros:
        1 zero acceleration/ 2 zeros velocity/ 3 zeros displacement

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband translation.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rt for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of ac for p-coda calculation.
    :type station: str
    :param station: Station of interest.
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param startev: Event origin time
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Detrended, trimmed rotation signal.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Detrended, trimmed three component translation.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: Detrended and trimmed copy rt for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: Detrended and trimmed copy of ac for p-coda calculation.
    """

    # translation poles and zeros dictionaries, output units of nm/s^2
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

        # different translation instruments require different method
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
        raise RotationalProcessingException("Invalid station")

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
    g = (np.sqrt(beta/np.pi) * 
         np.exp(-beta * (np.abs(freq - freq0) / bandwidth) ** 2.)) 

    # convolve your signal by the filter in frequency domain 
    sigarray_fourier = fft(sigarray) 
    sigarray_fourier_filtered = sigarray_fourier * g

    # back to time domain
    sigarray_filtered = np.real(ifft(sigarray_fourier_filtered))
    sigarray_filtered = highpass(sigarray_filtered, freq=0.0033, 
                                                df=5, corners=3, zerophase=True)
    sigarray_filtered = detrend(sigarray_filtered)

    return sigarray_filtered


def filter_and_rotate(rt, ac, rt_pcoda, ac_pcoda, cutoff, cutoff_pc, is_local):

    """
    Filters streams: lowpass at cutoff, and highpass with zerophase filter.
    Creates copies of ac/rt for pcoda analysis, both low and high sampling rate.
    Create lists of stream copies for analysis in different frequency bands.
    Output all new streams for use in later processing and plotting.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband translation.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rt for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of ac for p-coda calculation.
    :type cutoff/cutoff_pc: float
    :param cutoff/cutoff_pc: Cut-off frequency for lowpass filter.
    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :rtype trv_acc: :class: `~obspy.core.stream.Stream`
    :return trv_acc: T comp. of copied ac, low/highpass and bandpassed
    :rtype trv_pcoda: :class: `~obspy.core.stream.Stream`
    :return trv_pcoda: Copy of T component of ac_pcoda, not filtered
    :rtype rt_bands: list of streams, :class: `~obspy.core.stream.Stream`
    :return rt_bands: Copied rt, bandpass filtered at various freqs.
    :rtype trv_bands: list of streams, :class: `~obspy.core.stream.Stream`
    :return trv_bands: T comp. of copied ac, bandpass filtered at various freqs.
    :rtype rt_pcoda_coarse: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda_coarse: Copy of rt, highpass filtered at cutoff_pc
    :rtype trv_pcoda_coarse: :class: `~obspy.core.stream.Stream`
    :return trv_pcoda_coarse: T comp. of ac, highpass filtered at cutoff_pc
    :rtype filt_rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return filt_rt_pcoda: Copy of rt_pcoda, highpass filtered at cutoff_pc
    :rtype filt_ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return filt_ac_pcoda: Copy of ac_pcoda, highpass filtered at cutoff_pc
    :rtype filt_trv_pcoda: :class: `~obspy.core.stream.Stream`
    :return filt_trv_pcoda: T comp. of ac_pcoda, highpass filtered at cutoff_pc
    """

    # set the list of frequencies for bandpass filters
    freq_list = [0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0]
    number_of_bands = len(freq_list) - 1
    
    # lower sampling rate copies for pcoda analysis in page 2
    ac_pcoda_coarse = ac.copy()
    rt_pcoda_coarse = rt.copy()

    # filter base streams high and low
    ac.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    rt.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    ac.filter('highpass', freq=0.005, corners=2, zerophase=True)
    rt.filter('highpass', freq=0.005, corners=2, zerophase=True)
    ac_pcoda_coarse.filter(
                        'highpass', freq=cutoff_pc, corners=2, zerophase=True)
    rt_pcoda_coarse.filter(
                        'highpass', freq=cutoff_pc, corners=2, zerophase=True)

    # filter out secondary microseisms (5-12s) for far events
    if is_local == "far":    
        ac.filter('bandstop', freqmin=1/12, freqmax=1/5, 
                                corners=4, zerophase=True)
        rt.filter('bandstop', freqmin=1/12, freqmax=1/5, 
                                corners=4, zerophase=True)

    # single out transverse acceleration for processing
    trv_acc_tmp = ac.copy()
    trv_acc = trv_acc_tmp.rotate(method = 'NE->RT').select(component='T')

    # rotate pcoda streams to theoretical event backazimuth, for use in page 4
    pcoda_rotate = ac_pcoda.copy()
    pcoda_rotate_coarse_tmp = ac_pcoda_coarse.copy()

    trv_pcoda = pcoda_rotate.rotate(method = 'NE->RT').select(component='T')
    trv_pcoda_coarse = pcoda_rotate_coarse_tmp.rotate(method = 'NE->RT').select(
                                                                component='T')

    # highpass filter pcoda (w/ higher sampling rate), for use in page 4
    filt_rt_pcoda = rt_pcoda.copy()
    filt_ac_pcoda = ac_pcoda.copy()
    filt_ac_pcoda.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    filt_rt_pcoda.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)

    filt_ac_pcoda_tmp = filt_ac_pcoda.copy()
    filt_trv_pcoda = filt_ac_pcoda_tmp.rotate(method = 'NE->RT').select(
                                                                component='T')
    
    # for phase velocity estimation of varying frequency bands
    # rotate to TRansVerse component, filter both translations and rotations
    rt_bands = [rt.copy() for _ in range(number_of_bands)]
    trv_bands = [ac.copy() for _ in range(number_of_bands)]
    for I in range(number_of_bands):
        trv_bands[I] = trv_bands[I].rotate(method='NE->RT').select(
                                                                component='T')
        for band_select in [rt_bands[I],trv_bands[I]]:
            band_select[0].filter(type = 'bandpass',
                                freqmin = freq_list[I],
                                freqmax = freq_list[I+1],
                                corners = 3,
                                zerophase = True)

    return trv_acc, trv_pcoda, rt_bands, trv_bands, rt_pcoda_coarse, \
                trv_pcoda_coarse, filt_rt_pcoda, filt_ac_pcoda, filt_trv_pcoda 



def ps_arrival_times(ds_in_km, depth, init_sec):

    """
    Obtains the arrival times (in seconds after the start time of the fetched
    data) of the first P an S waves of the event. The inputs are the
    epicentral distance in degrees, the depth in km and the initial time in
    seconds (starttime_of_the_event - data_starttime)

    :type ds_in_km: float
    :param ds_in_km: Event-station distance, in km.
    :type depth: float
    :param depth: Hypocenter depth in km.
    :type init_sec: float
    :param init_sec: Initial time of the event in seconds.
    :rtype arriv_p: float
    :return arriv_p: P-wave first arrival.
    :rtype arriv_s: float
    :return arriv_s: S-wave first arrival.
    """
    # use taup to get the theoretical arrival times for P & S
    TauPy_model = TauPyModel('iasp91')
    tt = TauPy_model.get_travel_times(
                                distance_in_degree=ds_in_km / 111.11, 
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


def time_windows(dist_in_km, arriv_p, arriv_s, init_sec, is_local):

    """
    Determines time windows for P-waves S-waves, initial, latter surface waves. 
    Window lengths depends on event distance.
    All values set to integers, as they are used for slice indexing

    :type dist_in_km: float
    :param dist_in_km: Event station distance in km
    :type arriv_p: float
    :param arriv_p: Arrival time of the first P-wave.
    :type arriv_s: float
    :param arriv_s: Arrival time of the first S-wave.
    :type init_sec: float
    :param init_sec: Initial time of the event in sec in the fetched data.
    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :rtype `all`: int
    :return min_pw: Start time for P-waves window.
    :return max_pw: End time for P-waves window.
    :return min_sw: Start time for S-waves window.
    :return max_sw: End time for S-waves window.
    :return min_lwi: Start time for initial surface-waves window.
    :return max_lwi: End time for initial surface-waves window.
    :return min_lwf: Start time for latter surface-waves window.
    :return max_lwf: End time for latter surface-waves window.
    """

    if is_local == 'FAR':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + (arriv_s - arriv_p) // 4)
        min_sw = int(round(arriv_s - 0.001 * (arriv_s - arriv_p)))
        max_sw = int(arriv_s + 150)
        min_lwi = int(round(surf_tts(dist_in_km, init_sec) - 20))
        max_lwi = int(min_lwi + round((dist_in_km/1E3) * 50)) # 50sec/1000 km. 
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + round((dist_in_km/1E3) * 60)) # 60sec/1000 km.
    elif is_local == 'LOCAL':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + 20)
        min_sw = int(arriv_s - 5)
        max_sw = int(min_sw + 20)
        min_lwi = int(round(surf_tts(dist_in_km, init_sec) + 20))
        max_lwi = int(min_lwi + 50)
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + 80)
    elif is_local == 'CLOSE':
        min_pw = int(arriv_p)
        max_pw = int(min_pw + 7)
        min_sw = int(arriv_s)
        max_sw = int(min_sw + 7)
        min_lwi = int(round(surf_tts(dist_in_km, init_sec) + 5))
        max_lwi = int(min_lwi + 12)
        min_lwf = int(max_lwi)
        max_lwf = int(min_lwf + 80)


    return min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf


def surf_tts(ds_in_km, start_time):

    """
    Uses arrival times for different epicentral distances based on the IASP91
    travel times model to estimate a curve of travel times for surface waves.
    Returns the arrival time of the surface waves

    :type ds_in_km: float
    :param ds_in_km: Epicentral distance in km
    :type start_time: float
    :param start_time: Start time of the event.
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
        dife_r = abs(ds_in_km / 111.11 - np.arange(0., 180.1, 0.01)
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
    sw_arrival = arriv_lov + peq

    return sw_arrival


def get_corrcoefs(streamA, streamB, sec):

    """
    Calculates the zero-lag correlation coefficients between two streams in 
    small time windows, whose length is dictated by the 'sec' parameter
    *StreamA and StreamB data lengths need to be the same.

    :type streamA: :class: `~obspy.core.stream.Stream`
    :param streamA: First stream to correlate
    :type streamB: :class: `~obspy.core.stream.Stream`
    :param streamB: Second stream to correlate
    :type sec: int
    :param sec: Time window length.
    :rtype corrcoefs: numpy.ndarray
    :return corrcoefs: Correlation coefficients.
    :rtype thres: numpy.ndarray
    :return thres: Array for plotting dashed line of 75% correlation.
    """

    # time window in samples
    strA_TW = int(streamA[0].stats.sampling_rate * sec)
    strB_TW = int(streamB[0].stats.sampling_rate * sec)

    corrcoefs = []
    for i in range(0, len(streamA[0].data) // strA_TW):
        coeffs = correlate(a = streamA[0].data[i*strA_TW:(i+1)*strA_TW],
                           b = streamB[0].data[i*strB_TW:(i+1)*strB_TW], 
                           shift = 0)

        corrcoefs.append(coeffs[0])

    corrcoefs = np.asarray(corrcoefs)
    thres = 0.75 * np.ones(len(corrcoefs) + 1)

    return corrcoefs, thres


def baz_analysis(rt, ac, sec):

    """
    Computes correlation coefficients for varying backazimuth steps.
    Loops over BAz values from 0 to 360 in steps, for each BAz value, calculates
    correlations of rotation rate and tranvserse acceleration in small time 
    windows, whose length is controlled by the parameters 'sec'.
    *Data length of rt and ac need to be the same 

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Stream of the rotation data.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Stream of translation data.
    :type sec: int
    :param sec: Time window length.
    :rtype corrbaz_list: numpy.ndarray
    :return corrbaz_list: Array of correlation coefficients per backazimuth step
    :rtype maxcorr_list: numpy.ndarray
    :return maxcorr_list: BAz values for the maximum correlation per time window 
    :rtype backas: numpy.ndarray
    :return backas: Vector containing backazimuths values by step length
    :rtype coefs_list: list
    :return coefs_list: List of the max correlation value for each time window.
    """

    # time window in samples
    rt_TW = int(rt[0].stats.sampling_rate * sec)
    ac_TW = int(ac[0].stats.sampling_rate * sec)

    # separate data from streams for faster rotation/correlation
    acN = ac.select(component='N')[0].data
    acE = ac.select(component='E')[0].data
    rtZ = rt.select(component='Z')[0].data

    # create a list of backazimuths to iterate over
    step = 10
    corr_length = len(rt[0].data) // rt_TW
    backas = np.linspace(0, 360 - step, 360 / step)
    
    # iterate over BAz, rotate, correlate trace
    corrbaz_list = []
    for BAZ in backas:
        for j in range(0, corr_length):
            acRTZ = rotate_ne_rt(n = acN, e = acE, ba = BAZ)
            acT = acRTZ[1]
            corrbaz = correlate(a = rtZ[j*rt_TW:(j+1)*rt_TW],
                                b = acT[j*ac_TW:(j+1)*ac_TW],
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


def estimate_baz(rt, ac, start, end):

    """
    Estimate the backazimuth of an event by taking the average of all BAz's 
    which give a correlation greater than 0.9 in the s-wave and surface waves

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband translation
    :type start: float
    :param start: Starttime for S-waves window.
    :type end: float
    :param end: Endtime for latter surface waves window.
    :rtype corrsum_list: list of floats
    :return corrsum_list: Sum of all corr. coefficients above 0.9, start to end
    :rtype baz_list: numpy.ndarray
    :return baz_list: Vector containing backazimuths by step length
    :rtype max_ebaz_xcoefs: numpy.ndarray
    :return max_ebaz_xcoefs: Array of maximum correlations for each est. BAz
    :rtype EBA: float
    :return EBA: The estimated BAz if applicable, else NaN
    """

    # set integer sampling rates
    rt_SR = int(rt[0].stats.sampling_rate)
    ac_SR = int(ac[0].stats.sampling_rate)
    
    # time window in samples
    sec_internal = 30
    rt_TW = sec_internal * rt_SR
    ac_TW = sec_internal * ac_SR

    # sample number of surface wave start/end
    start_sample = start * rt_SR
    end_sample = end * ac_SR

    # cut streams at surface waves
    rt_cut = rt[0].data[start_sample:end_sample]
    acN_cut = ac.select(component='N')[0].data[start_sample:end_sample]
    acE_cut = ac.select(component='E')[0].data[start_sample:end_sample]


    step = 1
    baz_list = np.linspace(0, int(360 - step), int(360 / step)) # BAz array
    
    # iterate over all BAz values, correlate in time windows
    corr_list = []
    for BAZ in baz_list:
        for j in range(len(rt_cut) // rt_TW):
            rad_cut,trv_cut = rotate_ne_rt(n = acN_cut, e = acE_cut, ba = BAZ)
            corr = correlate(a = rt_cut[j*rt_TW:(j+1)*rt_TW],
                             b = trv_cut[j*ac_TW:(j+1)*ac_TW],
                             shift = 0)  

            corr_list.append(corr[0])

    corr_list = np.asarray(corr_list)
    corr_list = corr_list.reshape(
                        len(baz_list),int(round(len(corr_list)/len(baz_list))))

    # iterate over correlations, choose those > 0.9
    corrsum_list = []
    for k in range(len(baz_list)):
        bazsum = []
        for l in range(len(corr_list[0, :])):
            if corr_list[k, l] >= 0.9:
                bazsum.append(corr_list[k, l])
            else:
                bazsum.append(0.0)

        bazsum = np.asarray(bazsum)
        bazsum = sum(bazsum)
        corrsum_list.append(bazsum)

    # determine estimated backazimuth
    best_ebaz = baz_list[np.asarray(corrsum_list).argmax()] 
    max_ebaz_xcoef = np.max(corr_list[int(best_ebaz)]) 

    if max(corrsum_list) == 0.0:
        EBA = np.nan
    else:
        EBA = best_ebaz

    return corrsum_list, baz_list, max_ebaz_xcoef, EBA


def get_phase_vel(rt, trv_acc, sec, corrcoefs, start):

    """
    Calculate phase velocities by taking amplitude ratios, only for 
    correlation values > 0.75. 'start' controls where in data calculations begin

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Tranvserse acceleration stream.
    :type sec: int
    :param sec: Time window length.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Calculated correlation coefficients.
    :type start: int
    :param start: index from which to start calculations
    :rtype phasv_list: numpy.ndarray
    :return phasv_list: Calculated phase velocities
    """
    # time window in samples
    rt_TW = int(rt[0].stats.sampling_rate * sec)
    trv_TW = int(trv_acc[0].stats.sampling_rate * sec)

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


def sn_ratio(stream, p_arrival):

    """
    Characterizes the signal-to-noise ratio as the ratio of max amplitude and 
    the average signal value in the noise window before the first theoretical 
    p-wave arrival, assuming that the noise has the same behavior everywhere. 

    :type full_signal: numpy.ndarray
    :param full_signal: Amplitude data of the full signal.
    :type p_arrival: float
    :param p_arrival: Arrival time of the first P-wave.
    :type sam_rate: float
    :param sam_rate: Sampling rate.
    :rtype SNR: float
    :return SNR: Signal-to-noise ratio of the seismogram.
    """
    # convert to integers for indexing
    str_SR = int(stream[0].stats.sampling_rate)
    p_arrival = int(round(p_arrival))

    # determine signal max and noise value
    noise_data = stream[0].data[(p_arrival-180)*str_SR:(p_arrival-100)*str_SR]
    noise_mean = abs(np.mean(noise_data))
    data_max = max(stream[0].data)

    SNR = data_max/noise_mean

    return SNR


def store_info_json(rt, ac, trv_acc, data_sources, station, event, dist_baz, 
                    arriv_p, corrcoefs, EBA, max_ebaz_xcoef, 
                    phasv_means, phasv_stds, folder_name, tag_name):

    """
    Generates a human readable .json file to store processed data for each event

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband translations
    :type trv_acc: :class: `~obspy.core.stream.Stream`
    :param trv_acc: Transverse acceleration stream
    :type data_sources: dictionary
    :param data_sources: collection of data source for each channel.
    :type station: str
    :param station: Station of interest.
    :type event: :class: `~obspy.core.event.event.Event`
    :param event: Event information container
    :type dist_baz: tuple
    :param dist_baz: [0] Great circle distance in m, 
                [1] azimuth A->B in degrees,
                [2] backazimuth B->A in degrees.
    :type arriv_p: float
    :param arriv_p: P-wave first arrival time.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type EBA: float
    :param EBA: Estimated backazimuth.
    :type max_ebaz_xcoef: numpy.ndarray
    :return max_ebaz_xcoefs: Array of maximum correlations for each est. BAz
    :type phasv_means: numpy.ndarray
    :param phasv_means: Vector of mean phase velocities per freq. band
    :type phasv_std: numpy.ndarray
    :param phasv_std: Vector of phase velocities std. per freq. band
    :type folder_name: string
    :param folder_name: Name of the folder containing the event.
    :type tag_name: string
    :param tag_name: Handle of the event.
    """
    # set the global 'round to decimal point' value
    rnd = 6

    # parse out parameters for json file
    orig = event.preferred_origin() or event.origins[0] # Event origin
    catalog = orig.creation_info.author or orig.creation_info.agency_id
    magnitude = event.preferred_magnitude() or event.magnitudes[0] # Mag info.

    PAT = round(max(trv_acc[0].data), rnd)  # Peak transverse acc. [nm/s]
    PRZ = round(max(rt[0].data), rnd)  # Peak vertical rotation rate [nrad/s]
    PCC = round(max(corrcoefs), rnd)  # Peak correlation coefficient
    MCC = round(min(corrcoefs), rnd) # Minimum correlation coefficient
    TBA = round(dist_baz[2], rnd) # Theoretical backazimuth [°]
    MXE = round(max_ebaz_xcoef, rnd) # Max correlation for Estimated BAz
    DS_KM = round(0.001 * dist_baz[0], rnd) # Epicentral Distance [km]
    DS_DEG = round(DS_KM / 111.11, rnd) # Epicentral Distance [°]
    DS_CAT = is_local(DS_KM).lower()
    SNT = round(sn_ratio(ac, arriv_p), rnd)
    SNR = round(sn_ratio(rt, arriv_p), rnd)

    phasv_means = [round(_,rnd) for _ in phasv_means] 
    phasv_stds = [round(_,rnd) for _ in phasv_stds] 

    # common event dictionary
    dic_event = OrderedDict([
                ('event_id', event.resource_id.id),
                ('event_source', catalog),
                ('event_latitude', orig.latitude),
                ('event_longitude', orig.longitude),
                ('origin_time', str(orig.time)),
                ('trace_start', str(orig.time-180)),
                ('trace_end', str(orig.time+3*3600)),
                ('magnitude', magnitude.mag),
                ('magnitude_type', magnitude.magnitude_type),
                ('depth', orig.depth * 0.001),
                ('depth_unit', 'km')
                ])

    # individual station dictionary w/ rotational parameters and velocities
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
                    ('distance_category', DS_CAT),
                    ('theoretical_backazimuth', TBA),
                    ('estimated_backazimuth', EBA),
                    ('backazimuth_unit', 'degrees'),
                    ('peak_correlation_coefficient', PCC),
                    ('minimum_correlation_coefficient', MCC),
                    ('max_xcoef_for_estimated_backazimuth', MXE),
                    ('peak_vertical_rotation_rate', PRZ),
                    ('peak_vertical_rotation_rate_unit', 'nrad/s'),
                    ('peak_transverse_acceleration', PAT),
                    ('peak_transverse_acceleration_unit', 'nm/s^2'),
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
    json.dump(dic_event, outfile, indent = 4)
    outfile.close()


def store_info_xml(event,folder_name,tag_name,station):

    """
    Write QuakeML file. Store extra parameters under the namespace rotational
    seismology. Stations are taken care of in nested tags in the extra tag
    of the XML file. Parameters used for filtering events on JANE database 
    framework are taken from the previously written .json file 

    :type event: :class: `~obspy.core.event.event.Event`
    :param event: Event information container
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


def plot_waveform_comp(event, station, mode, folder_name, tag_name):

    """
    Main processing script, calls all other functions defined above.
    Compare vertical rotation rate and transverse acceleration.
    Creates and saves four figures, a .json file with processed parameters,
    and a QuakeML file for each event.
 
    :type event: :class: `~obspy.core.event.event.Event`
    :param event: Event information container
    :type station: str
    :param station: Station from which data are fetched (i.e. 'RLAS').
    :type link: string
    :param link: URL to the IRIS-XML file
    :type mode: str
    :param mode: Determines if moment tensor information is fetched
    :type folder_name: string
    :param folder_name: Name of the folder containing the event.
    :type tag_name: string
    :param tag_name: Handle of the event.
    """
    # =========================================================================
    #                                   
    #              Gather event information, stream objects etc.
    #
    # =========================================================================
    event_lat, event_lon, depth, startev, rt, ac, dist_baz, data_sources = \
        event_info_data(event, station, polarity, instrument)
    
    # parse out event and station location information
    ds_in_km = dist_baz[0] * 1E-3
    BAz = dist_baz[2]
    station_lat = rt[0].stats.coordinates.latitude
    station_lon = rt[0].stats.coordinates.longitude

    # parse out information for title
    flinn_engdahl_title = event.event_descriptions[0]['text'].title()
    mag = event.preferred_magnitude() or event.magnitudes[0] 

    # try to get moment tensor information
    if event.preferred_focal_mechanism():
        moment_tensor = get_moment_tensor(event)
    else:
        moment_tensor = None

    # =========================================================================
    #                                   
    #                                  PAGE 1
    #               Create map with event location & information
    #
    # =========================================================================
    print("\nPage 1 > Title Card...", end=" ")

    # ================================ Draw Maps===============================
    if is_local(ds_in_km) == 'FAR': 
        shift_text = 200000
        event_size = 200
        if ds_in_km <= 13000:
            plt.figure(figsize=(18, 9))
            plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

            # globe plot
            map = Basemap(projection='ortho', 
                            lat_0=(station_lat + event_lat) / 2, 
                            lon_0=(station_lon + event_lon) / 2, 
                            resolution='l')
            map.drawmeridians(np.arange(0, 360, 30))
            map.drawparallels(np.arange(-90, 90, 30))
        elif ds_in_km > 13000:
            plt.figure(figsize=(18, 9))
            plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

            # If the great circle between the station and event is crossing the
            # 180° meridian in the pacific and the stations are far apart the 
            # map has to be re-centered, otherwise wrong side of globe.
            if abs(station_lon - event_lon) > 180:
                lon0 = 180 + (station_lon + event_lon) / 2
            else:
                lon0 = (station_lon + event_lon) / 2

            map = Basemap(projection='moll', lon_0=lon0, resolution='l')
            map.drawparallels(np.arange(-90, 90, 30))
            map.drawmeridians(np.arange(0, 360, 30))

    elif is_local(ds_in_km) == 'LOCAL':
        shift_text = 35000
        event_size = 300
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

    elif is_local(ds_in_km) == 'CLOSE':
        event_size = 300
        plt.figure(figsize=(26, 13))
        plt.title('{}T{}Z\n \n '.format(startev.date, startev.time), 
                                                fontsize=24, fontweight='bold')
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

    if is_local(ds_in_km) == 'LOCAL' or is_local(ds_in_km) == 'CLOSE':
        map.drawlsmask(land_color='coral', ocean_color='lightblue', lakes=True)
        map.drawcountries(linewidth=0.6)
   
    map.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, 
                                                    linewidth=3, color='yellow')
    
    # =========================== Station/ Event ===============================
    ev_x, ev_y = map(event_lon, event_lat)
    sta_x, sta_y = map(station_lon, station_lat)

    # station
    map.scatter(sta_x, sta_y, 200, color='b', marker='v',
                                    edgecolor='k', zorder=100)
    plt.text(sta_x + shift_text, sta_y, station, va='top',
                                             family='monospace', 
                                             weight='bold', 
                                             zorder=101,
                                             color='k', 
                                             backgroundcolor='white')
    # event as moment tensor or star 
    # !!! doesn't work for some reason - disregard for now
    # if moment_tensor:
    #     ax = plt.gca()
    #     b = beach(moment_tensor, xy=(ev_x,ev_y), facecolor='blue',
    #                                 width=100, linewidth=1, alpha=1.0)
    #     b.set_zorder(100)
    #     ax.add_collection(b)
    # else:
    #     map.scatter(ev_x, ev_y, 200, color="b", marker="*", 
    #                                         edgecolor="k", zorder=200)

    # plot event
    map.scatter(ev_x, ev_y, event_size, color="b", marker="*", 
                                            edgecolor="k", zorder=100)

    # title large
    plt.subplot2grid((4, 9), (1, 0), colspan=2)
    plt.title(u'{}T{}Z\n'.format(startev.date, startev.time), 
                                        fontsize=20, weight='bold')
    ax = plt.gca()
    ax.axis('equal')
    ax.axis('off')

    # sub-title
    plt.subplot2grid((4, 9), (2, 0), colspan=2)
    plt.title(u'\n\nRegion: {}'.format(flinn_engdahl_title) + 
            '\n\nMagnitude: {} {}'.format(mag.mag,mag.magnitude_type) + 
            '\n\nDistance: {} [km], {} [°]'.format(
                                round(ds_in_km,2), round(BAz,2)) + 
            '\n\nDepth: {} [km]'.format(depth),
              fontsize=18, fontweight='bold')

    ax = plt.gca()
    ax.axis('off')

    plt.subplot2grid((4, 9), (3, 0), colspan=2)
    plt.title(u'Event Information: \n Global Centroid-Moment-Tensor '
              'Catalog (GCMT) \n\n Processing Date:\n' + 
              str(UTCDateTime().date),
              fontsize=14)

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
                                                    is_local(ds_in_km), rt, ac)

    print("Removing instrument response...")
    # remove instrument response based on station
    rt, ac, rt_pcoda, ac_pcoda = remove_instr_resp(
                                    rt, ac, rt_pcoda, ac_pcoda,station, startev)

    print("Filtering and rotating traces...")
    # filter raw data, rotate some to theoretical backazimuth, separate Pcoda
    trv_acc, trv_pcoda, rt_bands, trv_bands, rt_pcoda_coarse, trv_pcoda_coarse,\
    filt_rt_pcoda, filt_ac_pcoda, filt_trv_pcoda = filter_and_rotate(
            rt, ac, rt_pcoda, ac_pcoda, cutoff, cutoff_pc, is_local(ds_in_km))

    print("Getting theoretical arrival times...")
    # find trace start
    init_sec = startev - ac[0].stats.starttime

    # theoretical arrival times for P and S waves
    arriv_p, arriv_s = ps_arrival_times(ds_in_km, depth, init_sec)
    
    # determine time windows for seismic phases (P,S,surface)
    min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf = \
        time_windows(ds_in_km, arriv_p, arriv_s, init_sec, is_local(ds_in_km))

    # ======================================================================== 
    #                                
    #                                   Page 2
    #           Waveform Comparison plot w/ individual phase subplots
    #
    # =========================================================================    
    print("\nPage 2 >  Waveform Comparison...",end=" ")

    # rt.taper(max_percentage=0.05) # this was here but we already taper?

    # phase velocity, factor for displacing rotation rate, and time array
    c1 = .5 * max(abs(trv_acc[0].data)) / max(abs(rt[0].data))  
    fact1 = 2 * max(rt[0].data)  
    time = rt[0].stats.delta * np.arange(0, len(rt[0].data))

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
    
    # place box in middle of figure
    box_yposition = ((fact1 + max((1. / (2. * c1)) * trv_acc[0].data))
                     - abs(min(rt[0].data)))/2  
    
    # gap between annotation and vertical
    if is_local(ds_in_km) == 'FAR':  
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
    
    # pwave arrival times in samples for normal trace and pcoda
    min_pw_rt = rt_SR * min_pw
    max_pw_rt = rt_SR * max_pw
    min_pw_pcrt = rt_pcoda_SR * min_pw
    max_pw_pcrt = rt_pcoda_SR * max_pw

    # pcoda phase velocity
    cp = 0.5 * (max(abs(trv_pcoda_coarse[0][min_pw_pcrt:max_pw_pcrt]))/
                max(abs(rt_pcoda_coarse[0].data[min_pw_pcrt:max_pw_pcrt])))

    # find min and max trace amplitudes for y-limits
    min_ta_pcod = min((0.5 / cp) * trv_pcoda_coarse[0][min_pw_rt:max_pw_rt])
    max_ta_pcod = max((0.5 / cp) * trv_pcoda_coarse[0][min_pw_rt:max_pw_rt])
    min_rt_pcod = min(rt_pcoda_coarse[0].data[min_pw_rt:max_pw_rt])
    max_rt_pcod = max(rt_pcoda_coarse[0].data[min_pw_rt:max_pw_rt])

    plt.plot(time, rt_pcoda_coarse[0].data, color='r')
    plt.plot(time, (0.5 / cp) * trv_pcoda_coarse[0], color='k')
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
    corrcoefs, thres = get_corrcoefs(rt, trv_acc, sec)

    # calculate correlations for different frequency bands,
    # length of time windows given by seconds_list
    corrcoefs_bands, thresholds = [], []
    seconds_list = [200, 100, 50, 20, 12, 10, 8, 6]

    for i in range(len(rt_bands)):
        corrcoefs_tmp, thresh_tmp = get_corrcoefs(streamA = rt_bands[i],
                                                  streamB = trv_bands[i],
                                                  sec = seconds_list[i])
        corrcoefs_bands.append(corrcoefs_tmp)
        thresholds.append(thresh_tmp)

    # zero-lag correlation coefficients for range of backazimuths
    print("Analyzing correlation by BAz bins...")
    corrbaz, maxcorr, backas, max_coefs_10deg = baz_analysis(rt, ac, sec)

    # estimate backazimuth and correlations for given BAz
    print("Estimating best backazimuth values...")
    corrsum, baz_list, max_ebaz_xcoef, EBA = estimate_baz(
                                                        rt, ac, min_sw, max_lwf)

    print("Calculating phase velocities...")
    phasv = get_phase_vel(rt, trv_acc, sec, corrcoefs, start=0)
    
    # calculate phase velocities for different frequency bands
    surf_start = min_lwi // sec
    phasv_bands,phasv_means,phasv_stds = [],[],[]
    for i in range(len(rt_bands)):
        phasv_tmp = get_phase_vel(rt_bands[i], trv_bands[i], seconds_list[i],
                                        corrcoefs_bands[i], start=surf_start)
        
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

    # ======================================================================== 
    #                                
    #                                Page 3
    #              Cross Correlation, phase velocity determination,
    #                    Estimation of backazimuth figures
    #
    # ========================================================================= 
    print("\nPage 3 > Cross-correlation, Phase velocity...",end=" ")

    X, Y = np.meshgrid(np.arange(0, sec * len(corrcoefs), sec), backas)

    # subplot 1
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

    # subplot 2
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
                            ticks=np.linspace(0.75,1,6).tolist())
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')

    # subplot 3
    plt.subplot2grid((4, 26), (2, 0), colspan=25)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), max_coefs_10deg, 'ro-', 
                                    label='Max. CC for est. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), corrcoefs, 'ko-', 
                                        label='CC for theo. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), thres, '--r', lw=2)
    plt.ylabel(r'X-corr. coeff.', fontsize=10, fontweight='bold')

    # 75% annotation location
    if is_local(ds_in_km) == 'FAR':
        shift75 = 50
    else:
        shift75 = 10

    plt.text(time[len(time) - 1] + shift75, 0.71, r'0.75', color='red')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))

    min_corr = min(min(max_coefs_10deg), min(corrcoefs))
    plt.ylim(min_corr, 1)
    plt.legend(loc=4, shadow=True)
    plt.grid(True)

    # subplot 4
    plt.subplot2grid((4, 26), (3, 0), colspan=25)
    teobaz = BAz * np.ones(len(corrcoefs) + 1)
    plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r, vmin=-1, vmax=1)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), teobaz, '--r', lw=2)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), maxcorr, '.k')
    plt.text(1000, BAz, str(BAz)[0:5] + r'°',
             bbox={'facecolor': 'black', 'alpha': 0.8}, color='r')

    # only plot estimated backazimuth if a value is given
    if not np.isnan(EBA):
        obsbaz = EBA * np.ones(len(corrcoefs) + 1)
        plt.text(400, EBA, str(EBA)[0:5] + r'°',
                 bbox={'facecolor': 'black', 'alpha': 0.8}, color='y')
        plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec),
                 obsbaz, '--y', lw=2)

    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.xlabel(r'Time [s]', fontweight='bold')
    plt.ylabel(r'BAz [°]', fontsize=10, fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks(np.linspace(0,360,7).tolist())
    plt.grid(True)
    fig = plt.subplot2grid((4, 26), (3, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(fig, cmap=plt.cm.RdYlGn_r, norm=norm, 
                                                        orientation='vertical')
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks(np.linspace(-1,1,9).tolist())

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
    print("Analyzing rotations in the P-coda...")
    
    # Zero-lag correlation coefficients
    print("Zero-lag cross correlations...",end=" ")

    # integer sampling rates
    rt_pc_SR = int(filt_rt_pcoda[0].stats.sampling_rate)
    ac_pc_SR = int(filt_trv_pcoda[0].stats.sampling_rate)
    
    corrcoefs_p = []
    lwi_average = int(round((min_lwi+max_lwi)/2))

    # separate vertical components
    acZ_pcoda = ac_pcoda.select(component='Z')

    # cut pcoda at correct time window and taper cuts
    rt_pcoda_cut = filt_rt_pcoda.copy()
    ac_pcoda_cut = filt_ac_pcoda.copy()
    trv_pcoda_cut = filt_trv_pcoda.copy()

    rt_pcoda_cut[0].data = filt_rt_pcoda[0].data[0:lwi_average * rt_pc_SR]
    trv_pcoda_cut[0].data = trv_pcoda_cut[0].data[0:lwi_average * ac_pc_SR]
    for i in range(3):
        ac_pcoda_cut[i].data = filt_ac_pcoda[i].data[0:lwi_average * ac_pc_SR]
    
    for traces in [rt_pcoda,rt_pcoda_cut,ac_pcoda_cut,trv_pcoda_cut]:
        traces.taper(max_percentage=0.05)

    # find correlations
    corrcoefs_p, thres_p = get_corrcoefs(rt_pcoda_cut, trv_pcoda_cut, sec_p)

    print("Backzimuths...")
    # surface wave start sample
    max_lwi_ac = ac_pc_SR * max_lwi

    # analyze backazimuth
    corrbaz_p, maxcorr_p, backas_p, max_coefs_10deg_p = baz_analysis(
                                            rt_pcoda_cut, ac_pcoda_cut, sec_p)

    # set up arrays for plotting
    time_p = rt_pcoda_cut[0].stats.delta * np.arange(0, len(rt_pcoda[0].data))
    fact1_p = 2 * max(rt_pcoda[0].data[0:max_lwi_ac])

    c1_p = .5 * (max(abs(trv_pcoda[0].data[0:max_lwi_ac])) /
                 max(abs(rt_pcoda[0].data[0:max_lwi_ac])))

    # check for correlations >= 0.5
    maxcorr_p_list = []
    for m in range(0, len(maxcorr_p)):
        if np.max(corrbaz_p[:,m]) >= 0.5:
            maxcorr_p_list.append(maxcorr_p[m])
        else:
            maxcorr_p_list.append(0)

    # ======================================================================== 
    #                                
    #                                Page 4
    #               Cross correlations for P-Coda time window
    #
    # ========================================================================= 
    # P-coda limits
    xlim1 = rt[0].stats.delta * len(rt[0].data) 
    Xp, Yp = np.meshgrid(np.arange(0,sec_p * len(corrcoefs_p), sec_p), backas_p)

    print("\nPage 4 > P-Coda Waveform Comparison...",end=" ")

    # subplot 1
    plt.figure(figsize=(18, 9))
    plt.subplot2grid((5, 26), (0, 0), colspan=25)
    plt.plot(time_p, acZ_pcoda[0].data, color='g')
    plt.ylabel(r'a$_\mathbf{Z}$ [nm/s$^2$]', fontweight='bold', fontsize=11)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(min(acZ_pcoda[0].data[0:max_lwi_ac]),
             max(acZ_pcoda[0].data[0:max_lwi_ac]))
    plt.title(r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ and a$_\mathbf{T}$'
              'correlation in the P-coda in a %d seconds time window'
              ' (highpass, cutoff: 1 Hz). Event: %s %sZ' % 
                                            (sec_p, startev.date, startev.time))
    plt.axvline(x=min_pw, linewidth=1)
    plt.axvline(x=min_sw, linewidth=1)
    plt.grid(True)

    # subplot 2
    xlim2 = (min_lwi + max_lwi) // 2
    plt.subplot2grid((5, 26), (1, 0), colspan=25, rowspan=2)

    plt.plot(time_p, rt_pcoda[0].data, color='r', label=r'Rotation rate')
    plt.plot(time_p, (0.5 / c1_p) * trv_pcoda[0].data + fact1_p, color='k',
                                             label=r'Transversal acceleration')
    plt.ylabel(r'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] -'
                   'a$_\mathbf{T}$/2c [1/s]', fontweight='bold', fontsize=11)
    plt.xlim(0, xlim2)
    plt.ylim(min(rt_pcoda[0].data[0:max_lwi_ac]), 
        fact1_p + max((1. / (2. * c1_p)) * trv_pcoda[0].data[0:max_lwi_ac]))
    
    box_yposition2 = (fact1_p + max((1. / (2. * c1_p)) * 
                        trv_pcoda[0].data[0:max_lwi_ac]) -
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

    # subplot 3
    plt.subplot2grid((5, 26), (3, 0), colspan=25)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p), corrcoefs_p, '.k')
    plt.ylabel(r'X-corr. coeff.', fontweight='bold')
    plt.xlim(0, xlim2)
    plt.ylim(0, 1)
    plt.grid(True)

    # subplot 4
    plt.subplot2grid((5, 26), (4, 0), colspan=25)
    plt.pcolor(Xp, Yp, corrbaz_p, cmap=plt.cm.RdYlGn_r, vmin=-1, vmax=1)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p), 
                                                        maxcorr_p_list, '.k')
    plt.xlim(0, xlim2)
    plt.xlabel(r'Time [s]', fontweight='bold')
    plt.ylabel(r'BAz [°]', fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks(np.linspace(0,360,7).tolist())
    plt.grid(True)

    fig = plt.subplot2grid((5, 26), (4, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(fig, cmap=plt.cm.RdYlGn_r, norm=norm, 
                                                        orientation='vertical')
    cb1.set_label(r'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks(np.linspace(-1,1,9).tolist())
   
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
    
    store_info_json(rt, ac, trv_acc, data_sources, station, event, dist_baz, 
                    arriv_p, corrcoefs, EBA, max_ebaz_xcoef, phasv_means, 
                    phasv_stds, folder_name, tag_name)

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
    :rtype check_folder_exists: list of str
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
        GCMT catalog for the most up to date catalog. ISC QuakeML file for \
        catalog of local/regional events. IRIS for most stable solutions, \
        though recent events might not be present \
        (default: gcmt, else: iscquakeml, iris)', type=str,default='GCMT')
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
        (default is 0 km). Positive down for IRIS.', 
                                            type=float or int, default=0.0)
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

    # [default]: get event catalog from GCMT NEW QUICK,
    if mode == 'GCMT':
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
        sys.exit('Invalid mode: {}\nValid: GCMT,IRIS,ISCQUAKEML'.format(mode))

    # set file output path
    output_path = './OUTPUT/'
    if not os.path.exists(output_path): 
        os.makedirs(output_path)

    print("%i event(s) downloaded, beginning processing...\n" % len(cat))
    event_counter = success_counter = fail_counter = already_processed = 0
    bars = '='*79
    error_list,error_type = [],[]
    for event in cat:
        event_counter += 1
        print("{} of {} event(s)".format(event_counter,len(cat)))
        try:
            tag_name, folder_name, check_folder_exists = generate_tags(event)
            # check if current event folder exists
            if check_folder_exists:
                # check if event source is the same, assumes 0 or 1 files found
                if (os.path.basename(check_folder_exists[0]) != 
                                                os.path.basename(folder_name)):
                    print("This event was processed with another mode\n")
                    error_list.append(tag_name)
                    error_type.append("Processed w/ Another Mode")
                    already_processed += 1
                    continue

                # if new station, run waveform compare again
                try:
                    filename_json = os.path.join(folder_name,tag_name + '.json')
                    data = json.load(open(filename_json))
                    if data['station_information_{}'.format(station)]:
                        print("This event was already processed\n")
                        error_list.append(tag_name)
                        error_type.append("Already Processed")
                        already_processed += 1
                    else:
                        try:
                            plot_waveform_comp(event, station, mode,
                                                        folder_name, tag_name)
                            success_counter += 1

                        # if any error, remove folder, continue
                        except Exception as e:
                            print(e)
                            print("Removing incomplete folder...\n")
                            error_list.append(tag_name)
                            error_type.append(e)
                            shutil.rmtree(folder_name)
                            fail_counter += 1

                        # if keyboard interrupt, remove folder, quit
                        except KeyboardInterrupt:
                            print("Removing incomplete folder...\n")
                            shutil.rmtree(folder_name)
                            sys.exit()

                # if json not found, folder is incomplete, continue
                except FileNotFoundError:
                    error_list.append(tag_name)
                    error_type.append("Incomplete Folder")
                    print("Incomplete folder found\n")
                    fail_counter += 1 


            
            # event encountered for the first time, create folder, xml, process
            elif not check_folder_exists:  
                os.makedirs(str(folder_name))
                
                # run processing function
                try:
                    plot_waveform_comp(event, station, mode, 
                                                        folder_name, tag_name)
                    success_counter += 1
                
                # if any error, remove folder, continue
                except Exception as e:
                    print(e)
                    print("Removing incomplete folder...\n")
                    error_list.append(tag_name)
                    error_type.append(e)
                    shutil.rmtree(folder_name)
                    fail_counter += 1

                # if keyboard interrupt, remove folder, quit
                except KeyboardInterrupt:
                    fail_counter += 1
                    print("Removing incomplete folder...\n")
                    shutil.rmtree(folder_name)
                    sys.exit()

        # if error creating tags, continue
        except Exception as e:
            print("Error in folder/tag name creation;\n",e)
            error_list.append(event.resource_id.id)
            error_type.append('Tag Creation')
            fail_counter += 1


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

        write_mode = 'w'
        if os.path.exists(errorlog_name):
            write_mode = 'a+'

        with open(errorlog_name, write_mode) as f:
            f.write("Error Log Created {}\n".format(datetime.datetime.now()))

            # prompt showing search parameters and number of failed processes
            f.write("{}<datetime<{}\n{}<mag<{}\nmode:{}\nfailed:{}/{}\n".format(
                                        args.min_datetime,args.max_datetime,
                                        args.min_magnitude,args.max_magnitude,
                                        mode,fail_counter+already_processed,
                                        len(cat)))
            for i,j in zip(error_list,error_type):
                f.write('{}\t{}\n'.format(i,j))
            f.write('_'*79)

# Debugger (* paste in wherever you want to break the code)
# import pdb; pdb.set_trace()
