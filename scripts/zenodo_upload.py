#!/usr/bin/env python
# some part are shamelessly stolen from https://github.com/darvasd/upload-to-zenodo

import os
import requests
import click
import json
import codecs

# global static header
headers = {"Content-Type": "application/json"}

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--filename',   '-f', type=str, help='Files to upload to Zenodo.org', required=False)
@click.option('--token',      '-t', type=str, help='Access token', required=True)
@click.option('--metadata',   '-m', type=str, help='Metadata file', required=False)
@click.option('--deposition', '-d', type=str, help='Deposition ID', required=False)
@click.option('--sandbox',    '-s', is_flag=True, help='Use sandbox', required=False, default=True)
# This is the surrogate for main(): everything happens here
def uploadToZenodo(filename,token,metadata,sandbox,deposition):
    print "Processing file ",filename
    zenodoHost='zenodo.org'
    if sandbox:
        zenodoHost='sandbox.' + zenodoHost
    print "Using host: " + zenodoHost

    # The URL for the actual upload and a global 
    zenodoURL='https://' + zenodoHost + '/api/deposit/depositions'

    # if there is no deposition ID, create one
    newDeposition = False
    if not deposition:
        deposition = new_deposition(zenodoURL, token)
        newDeposition = True
    # add file - in some cases we are only updating metadata, but 
    # to create a new deposition you must have a file
    if not newDeposition and filename is not None:
        upload_file(filename, deposition, zenodoURL, token)

    if metadata:
        print "Uploading metadata %s" % metadata
        r = upload_metadata(zenodoURL, metadata, deposition, token) 

def new_deposition(zURL,token):
    r = requests.post( zURL,
            params={'access_token': token}, 
            json={}, 
            headers=headers)
    depo_id = r.json()['id']
    print "New deposition created with ID %s at %s" % (depo_id, zURL +"/"+str(depo_id))
    print "Refer to this at further file uploads"
    return depo_id

def upload_file(filename, deposition_id, zURL, token):
    data = {'filename': os.path.basename(filename)}
    files = {'file': open(filename, 'rb')}
    r = requests.post(zURL + '/%s/files' % deposition_id,
            params={'access_token': token}, 
            data=data,
            files=files)
    # TODO: actualy it is not true: if the filename is already there, it is not uploaded, the JSON
    # returns with a warning and does nothing
    print "New file %s uploaded" % filename

def upload_metadata(zURL, mdFile, deposition_id, token):
     # reading metadata:
    with codecs.open(mdFile, 'r', 'utf-8') as f:
        metadataStr = f.read()
    if not _is_valid_json(metadataStr):
        return
    else:
        metadataJSON = json.loads(metadataStr)

    return requests.put(zURL + '/%s' % deposition_id, 
            params={'access_token': token}, 
            json=metadataJSON, 
            headers=headers)

def _is_valid_json(text):
    try:
        json.loads(text)
        return True
    except ValueError as e:
        print('Invalid json: %s' % e)
        return False

if __name__ == "__main__":
    uploadToZenodo()
