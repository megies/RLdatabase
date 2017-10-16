"""10.10.17
A script to interact with the database. Most useful for uploading single stationXML files
and deleting single events or stationxmls. Wipe is dangerous
"""
import os
import sys
import glob
import requests
import argparse

# command line arguments
parser = argparse.ArgumentParser(description='Database interaction script for'
    'adding stations, individual events, attaching individual attachments,'
    'deleting stations and events or choosing the nuclear option and wiping the'
    'database clean. All arguments needed except for wipe.',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--action', help='Performable actions:\n\
    get: return all XML information\n\
    put: upload a XML file\n\
    attach: post an attachment\n\
    delete: delete an event and or station\n\
    wipe: wipe all event XML files, no fileid', type=str, default='blank')
parser.add_argument('--pick',help='Choose [stationxml] or [quakeml]',type=str,
                                                                default='blank')
parser.add_argument('--fileid', help='Filename for object as input,\n' 
    'If [all], will post all stationxml files to database.\n'
    'For mass event upload, use event_upload_rotjane.py',
                                                    type=str,default='blank')
args = parser.parse_args()

action = args.action
pick = args.pick
filename = args.fileid
if (action or pick or filename) == 'blank':
    sys.exit('Requires all arguments, use -h to see options')


# ===============================MAIN==========================================
# jane settings
root_path = 'http://127.0.0.1:8000/rest/'
authority = ('chow','chow')

# folder with stataionXML files, quakeml files
sta_path = './station_files/'
eve_path = '/path/to/OUTPUT/'

# actions
if action == 'get':
    if filename == 'blank': filename = ''
    r = requests.get(
            url=root_path + 'documents/{}/{}'.format(pick,filename),
            auth=authority)

    assert r.ok
    print(r.content)

elif action == 'put':
    if filename == 'all' and pick == 'stationxml':
        # bit hacky because dataless in a folder, but works
        stations = glob.glob(sta_path + '*')
        station_names = [_[len(sta_path):] for _ in stations]
        for sta,sta_name in zip(stations,station_names):
            with open(sta,'rb') as fh:
                r = requests.put(
                url=root_path + 'documents/{}/{}'.format(pick,sta_name),
                auth=authority,
                data=fh)

    else:
        with open(filename,'rb') as fh:
            r = requests.put(
                url=root_path + 'documents/{}/{}'.format(pick,filename),
                auth=authority,
                data=fh)

    assert r.ok 


elif action == 'delete':
    r = requests.delete(
            url=root_path + 'documents/{}/{}'.format(pick,filename),
            auth=authority)

    assert r.ok

elif action == 'attach' and pick == 'quakeml':
    content = input('content-type?: ')    
    category = input('category?: ')
    headers = {'content-type': '{}'.format(content),
                 'category': '{}'.format(category)} 

    doc_ind = input('document index: ')
    Rtmp = requests.get(
            url=root_path + 'document_indices/{}/{}'.format(pick,doc_ind),
            auth=authority)
    print(Rtmp.content)

    with open(filename,'rb') as fh:
        r = requests.post(
            url=root_path + 'document_indices/{}/{}/attachments'.format(
                                                                pick,doc_ind),
            auth=authority,
            headers=headers,
            data=fh)

    assert r.ok

elif action == 'wipe':
    wipe = input('Are you sure you want to wipe the {} folder?\n'
                                                    'yes or no?: '.format(pick))
    if wipe == 'yes' and pick == 'quakeml':
        filelist = glob.glob(eve_path + 'GCMT*') + glob.glob(eve_path + 'ISC*')
        filenames = ['{}.xml'.format(_[len(eve_path):]) for _ in filelist]
        a=1/0
    elif wipe == 'yes' and pick == 'stationxml':
        filelist = glob.glob(sta_path + '*')
        filenames = [_[len(sta_path):] for _ in filelist]   
        a=1/0 
    else:
        sys.exit('Aborting')

    for filename in filelist:
        try:
            r = requests.delete(
                    url=root_path + 'documents/{}/{}'.format(pick,filename),
                    auth=authority)
            assert r.ok
        except AssertionError:
            continue
       



