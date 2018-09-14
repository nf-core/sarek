#!/usr/bin/env python
# some part are shamelessly stolen from https://github.com/darvasd/upload-to-zenodo

import os
import requests
import click
import json
import codecs

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--filename',   '-f', type=str, help='Files to upload to Zenodo.org', required=True)
@click.option('--token',      '-t', type=str, help='Access token', required=True)
@click.option('--metadata',   '-m', type=str, help='Metadata file', required=True)
@click.option('--deposition', '-d', type=str, help='Deposition ID', required=False)
@click.option('--sandbox',    '-s', is_flag=True, help='Use sandbox', required=False, default=True)

def uploadToZenodo(filename,token,metadata,sandbox,deposition):
    print "Processing file ",filename
    zenodoHost='zenodo.org'
    if sandbox:
        zenodoHost='sandbox.' + zenodoHost
    print "Using host: " + zenodoHost

    # reading metadata:
    with codecs.open(metadata, 'r', 'utf-8') as f:
        metadataStr = f.read()
    if not _is_valid_json(metadataStr):
        return
    else:
        metadataJSON = json.loads(metadataStr)
    print metadataJSON

    zenodoURL='https://' + zenodoHost + '/api/deposit/depositions'
    headers = {"Content-Type": "application/json"}
#    r = requests.post( zenodoURL,
#            params={'access_token': token}, 
#            json={}, 
#            headers=headers)

#    deposition_id = r.json()['id']
#    data = {'filename': os.path.basename(filename)}
#    files = {'file': open(filename, 'rb')}
#    r = requests.post(zenodoURL + '/%s/files' % deposition_id,
#            params={'access_token': token}, 
#            data=data,
#            files=files)

#    print r.json()

    deposition_id = deposition

    r = requests.put(zenodoURL + '/%s' % deposition_id, 
            params={'access_token': token}, 
            json=metadataJSON, 
            headers=headers)

    print r.json()


def _is_valid_json(text):
    try:
        json.loads(text)
        return True
    except ValueError as e:
        print('Invalid json: %s' % e)
        return False

if __name__ == "__main__":
    uploadToZenodo()
