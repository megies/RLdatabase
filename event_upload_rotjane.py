#!/usr/bin/env python3
"""04.10.17
For interacting with Rotational Jane database:
to upload events from event folders as quakeml files 
attach .png and .json files to each quakeml through REST format
"""
import re
import os
import sys
import glob
import argparse
import requests
import datetime
from pathlib import Path

# settings 
root_path = 'http://127.0.0.1:8000/rest/'
authority = ('chow','chow')
CURDIR = Path().absolute()
OUTPUT_PATH = CURDIR / 'OUTPUT'
ERRORLOGDIR = CURDIR / 'errorlogs'
# our apache is serving this via https now, so we have to use the Geophysik
# root certificate
# the certificate was switched to an official DFN certificate (that should be
# registered in the system)
# SSL_ROOT_CERTIFICATE = os.path.expanduser(
#     '~/ssl/geophysik_root_certificate/CAcert.pem')
requests_kwargs = {
    'auth': authority,
    # 'verify': SSL_ROOT_CERTIFICATE,
    }

# command line arguments
parser = argparse.ArgumentParser(description='Upload event quakeml and \
    post attachments to rotational Jane database.')
parser.add_argument('--timespan', help='What time span to upload files for, \
    options are the past [week] (default), [all] events available or instead \
    the path(s) to folder(s) for individual event(s) which overrides this \
    option.', type=str, default='week')
parser.add_argument(
    'directories', type=str, action="store", nargs="*", help='Individual \
    event folder(s) to upload. If specified, "timespan" mechanism is ignored.')
args = parser.parse_args()
timespan = args.timespan
directories = args.directories

# for attachments
png_re_pattern = re.compile(r'.*_([a-zA-Z]+)_page_([0-9]).png$')
png_page_category_map = {
    1: 'Event Information',
    2: 'Waveform Comparison',
    3: 'Correlation/Backazimuth',
    4: 'P-Coda Comparison'}
headers_json = {'content-type': 'text/json',
                'category': 'Processing Results'}

if directories:
    cat = list(Path(d).absolute() for d in directories)
    for path in cat:
        if not path.exists() or not path.is_dir():
            msg = (f'Location given on command line is not a directory: '
                   f'{str(path)}')
            raise ValueError(msg)
else:
    if timespan == 'week':
        # look for events in the past week
        cat = []
        for J in range(7):
            past = datetime.date.today() - datetime.timedelta(days=J)
            day = (OUTPUT_PATH / past.strftime('%Y') / past.strftime('%m')).glob(
                f'*_{past.isoformat()}_*')
            cat += list(day)
    elif timespan == 'all':
        # initial population, grab all events in folder
        cat = OUTPUT_PATH.glob('*/*/*_*_*')
    else:
        raise ValueError("bad 'timespan' option: '%s'" % timespan)
    
#cat = OUTPUT_PATH.glob('2024/06/*_*_*')

# ============================================================================

error_list,error_type = [],[]
for event in cat:
    os.chdir(CURDIR)
    print(event.name)
    try:
        os.chdir(event)
        attachments = sorted(glob.glob('*'))
        xml_files = [filename for filename in attachments
                     if filename.endswith('.xml')]

        if len(xml_files) != 1:
            error_list.append(event)
            error_type.append('No xml file for event')
            continue
        xml_file = xml_files[0]

        # check: full folder
        if len(attachments) < 6:
            error_list.append(event)
            error_type.append('Attachment Number Too Low: {}'.format(
                                                            len(attachments)))
            continue

        # push quakeml file
        with open(xml_file,'rb') as fh:
            r = requests.put(
                url=root_path + 'documents/quakeml/{}'.format(xml_file),
                data=fh, **requests_kwargs)

        # find attachment url
        r = requests.get(
                url=root_path + 'documents/quakeml/{}'.format(xml_file),
                **requests_kwargs)
        assert r.ok

        attachment_url = r.json()['indices'][0]['attachments_url']

        # delete all existing attachments, if there are any, before uploading
        # new attachments
        attachment_url_next = attachment_url
        r = requests.get(attachment_url_next, **requests_kwargs)
        assert r.ok
        while r.json()['count']:
            for result in r.json()['results']:
                print(f" deleting: {result['url']}")
                r_ = requests.delete(result['url'], **requests_kwargs)
                assert r_.ok
            r = requests.get(attachment_url, **requests_kwargs)
            assert r.ok

        # post image attachments            
        for filename in attachments:
            if filename.endswith('.xml'):
                continue
            print(filename)
            if filename.endswith('.png'):
                match = re.match(png_re_pattern, filename)
                if not match:
                    continue
                station = match.group(1)
                page_number = int(match.group(2))
                category = '{} ({})'.format(
                    png_page_category_map[page_number], station)
                if station == 'RLAS' or station == 'ROMY':
                    station = 'dataless.seed.BW_{}.xml'.format(station)
                else:
                    error_list.append(event)
                    error_type.append('Unidentified Station: {}'.format(station))
                    continue
                header = {'content-type': 'image/png',
                          'category': category,
                          'station': station}
            elif filename.endswith('.json'):
                header = headers_json
            else:
                error_list.append(event)
                error_type.append('Unidentified Attachment: {}'.format(J))

            with open(filename, 'rb') as fh:
                r = requests.post(url=attachment_url, headers=header, data=fh,
                                  **requests_kwargs)
            assert r.ok


        # the following is not used for now, the whole checking logic was hard
        # coded to a single station only and is pretty hacky.
        # XXX # check: already uploaded (409) and check for incomplete folders
        # XXX if r.status_code == 409:
        # XXX     r2 = requests.get(
        # XXX         url=root_path + 'documents/quakeml/{}'.format(xml),
        # XXX         **requests_kwargs)
        # XXX     assert r2.ok
        # XXX
        # XXX     try:
        # XXX         att_count = r2.json()['indices'][0]['attachments_count']
        # XXX         if att_count == 5:
        # XXX             continue
        # XXX         elif att_count != 5:
        # XXX             error_list.append(event)
        # XXX             error_type.append('Already Uploaded; Attachment Count Error')
        # XXX             continue
        # XXX     except IndexError:
        # XXX         error_list.append(event)
        # XXX         error_type.append('Already Uploaded; Attachment Count Error')
        # XXX         continue
        # XXX
        # XXX assert r.ok

    except ConnectionError:
        error_list.append(event)
        error_type.append('Connection Error')
        continue

    except AssertionError:
        # if assertion fails for any reason, delete current folder
        print(r.content.decode('UTF-8'))
        r_del = requests.delete(
                url=root_path + 'documents/quakeml/{}'.format(xml_file),
                **requests_kwargs)

        # tag errors for errolog
        error_list.append(os.path.basename(event))
        try:
            reason = r.json()['reason']
        except Exception:
            reason = ''
        error_type.append(str(r.status_code) + ' ' + reason)
        continue

# write error log to txt file to see what failed
if len(error_list) > 0:
    if not ERRORLOGDIR.exists():
        ERRORLOGDIR.mkdir()
    timestamp = datetime.datetime.now()
    M = timestamp.month
    Y = timestamp.year
    log_name = 'upload_{}_{}.txt'.format(M,Y)

    with open(ERRORLOGDIR / log_name, 'at') as f:
        f.write('Error Log Created {}\n'.format(timestamp))
        f.write('{}\n'.format('='*79)) 
        for i in range(len(error_list)):
           f.write('{}\n> {}\n'.format(error_list[i],error_type[i]))
        f.write('{}\n'.format('='*79)) 

    print('Logged {} error(s)'.format(len(error_list)))
