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

# settings 
root_path = 'http://127.0.0.1:8000/rest/'
authority = ('chow','chow')
OUTPUT_PATH = os.path.abspath('./OUTPUT/')
# our apache is serving this via https now, so we have to use the Geophysik
# root certificate
# the certificate was switched to an official DFN certificate (that should be
# registered in the system)
# SSL_ROOT_CERTIFICATE = os.path.expanduser(
#     '~/ssl/geophysik_root_certificate/CAcert.pem')
requests_kwargs = {
    'auth' : authority,
    # 'verify': SSL_ROOT_CERTIFICATE,
    }

# command line arguments
parser = argparse.ArgumentParser(description='Upload event quakeml and \
    post attachments to rotational Jane database.')
parser.add_argument('--timespan', help='What time span to upload files for, \
    options are the past [week] (default), or [all] events available.', 
                                                type=str, default='week')
args = parser.parse_args()
timespan = args.timespan

# for attachments
head_p1 = {'content-type': 'image/png',
                 'category': 'Event Information'} 
head_p2 = {'content-type': 'image/png',
                 'category': 'Waveform Comparison'} 
head_p3 = {'content-type': 'image/png',
                 'category': 'Correlation/Backazimuth'} 
head_p4 = {'content-type': 'image/png',
                 'category': 'P-Coda Comparison'} 
headers_json = {'content-type': 'text/json',
                 'category': 'Processing Results'} 

if timespan == 'week':
    # look for events in the past week
    cat = []
    for J in range(7):
        past = datetime.datetime.utcnow() - datetime.timedelta(days=J)
	#.  Look for folders inside two subdirectories      
        day = glob.glob(os.path.join(OUTPUT_PATH, '*/*/GCMT_{}*'.format(past.isoformat()[:10])))
        cat += day
elif timespan == 'all':
    # initial population, grab all events in folder
    cat = glob.glob(os.path.join(OUTPUT_PATH, '*/*/GCMT*')) + \
        glob.glob(os.path.join(OUTPUT_PATH, '*/*/ISC*'))
    cat.sort(reverse=True)
# ============================================================================

error_list,error_type = [],[]
for event in cat:
    print(os.path.basename(event))
    try:
        os.chdir(event)
        attachments = glob.glob('*')

        # check: full folder
        if len(attachments) < 6:
            error_list.append(event)
            error_type.append('Attachment Number Too Low: {}'.format(
                                                            len(attachments)))
            os.chdir('..')
            continue

        # assign attachments (kinda hacky)
        for J in attachments:
            if '.json' in J:
                json = J
            elif '.xml' in J:
                xml = J
            elif '_page_1.png' in J:
                page1 = J
            elif '_page_2.png' in J:
                page2 = J
            elif '_page_3.png' in J:
                page3 = J
            elif '_page_4.png' in J:
                page4 = J
            else:
                error_list.append(event)
                error_type.append('Unidentified Attachment: {}'.format(J))

        # push quakeml file
        with open(xml,'rb') as fh:
            r = requests.put(
                url=root_path + 'documents/quakeml/{}'.format(xml),
                data=fh, **requests_kwargs)

        # check: already uploaded (409) and check for incomplete folders
        if r.status_code == 409:
            r2 = requests.get(
                url=root_path + 'documents/quakeml/{}'.format(xml),
                **requests_kwargs)
            assert r2.ok

            try:
                att_count = r2.json()['indices'][0]['attachments_count']
                if att_count == 5:
                    os.chdir('..')
                    continue
                elif att_count != 5:
                    error_list.append(event)
                    error_type.append('Already Uploaded; Attachment Count Error')
                    os.chdir('..')
                    continue
            except IndexError:
                error_list.append(event)
                error_type.append('Already Uploaded; Attachment Count Error')
                os.chdir('..')
                continue

        assert r.ok

        # find attachment url
        r = requests.get(
                url=root_path + 'documents/quakeml/{}'.format(xml),
                **requests_kwargs)
        assert r.ok

        attachment_url = r.json()['indices'][0]['attachments_url']

        # post image attachments            
        for pngs,heads in zip([page1,page2,page3,page4],
                                [head_p1,head_p2,head_p3,head_p4]):
            with open(pngs,'rb') as fhp:
                r = requests.post(url=attachment_url, headers=heads, data=fhp,
                                  **requests_kwargs)

            assert r.ok

        # post .json
        with open(json,'rb') as fhj:
            r = requests.post(url=attachment_url, headers=headers_json,
                              data=fhj, **requests_kwargs)

            assert r.ok

        os.chdir('..')

    except ConnectionError:
        error_list.append(event)
        error_type.append('Connection Error')
        os.chdir('..')
        continue

    except AssertionError:
        # if assertion fails for any reason, delete current folder
        print(r.content.decode('UTF-8'))
        r_del = requests.delete(
                url=root_path + 'documents/quakeml/{}'.format(xml),
                **requests_kwargs)

        # tag errors for errolog
        error_list.append(os.path.basename(event))
        try:
            reason = r.json()['reason']
        except Exception:
            reason = ''
        error_type.append(str(r.status_code) + ' ' + reason)
        os.chdir('..')
        continue

# write error log to txt file to see what failed
if len(error_list) > 0:
    if not os.path.exists('../errorlogs'):
            os.makedirs('../errorlogs')
    os.chdir('../errorlogs')
    timestamp = datetime.datetime.now()
    M = timestamp.month
    Y = timestamp.year
    mode = 'w'
    log_name = 'upload_{}_{}.txt'.format(M,Y)

    # check if file exists
    if os.path.exists(log_name):
        mode = 'a+'

    with open(log_name,mode) as f:
        f.write('Error Log Created {}\n'.format(timestamp))
        f.write('{}\n'.format('='*79)) 
        for i in range(len(error_list)):
           f.write('{}\n> {}\n'.format(error_list[i],error_type[i]))
        f.write('{}\n'.format('='*79)) 


    print('Logged {} error(s)'.format(len(error_list)))


