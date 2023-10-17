import requests

headers = {"Content-Type": "application/json"}
params = {'access_token': ZENODO_DEPOSIT}

r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions',
                params=params,
                json={},
                '''
                Headers are not necessary here since "requests" automatically
                adds "Content-Type: application/json", because we're using
                the "json=" keyword argument
                headers=headers,
                '''
                headers=headers)
r.status_code
# 201
bucket_url = r.json()["links"]["bucket"]

''' New API '''
filename = "*.vcf.gz"
path = "./variant_calling/%s" % filename

'''
The target URL is a combination of the bucket link with the desired filename
seperated by a slash.
'''
with open(path, "rb") as fp:
    r = requests.put(
        "%s/%s" % (bucket_url, filename),
        data=fp,
        params=params,
    )
r.json()

data = {
    'metadata': {
        'title': 'My first upload',
        'upload_type': 'poster',
        'description': 'This is my first upload',
        'creators': [{'name': 'Doe, John',
                    'affiliation': 'Zenodo'}]
    }
}
r = requests.put('https://sandbox.zenodo.org/api/deposit/depositions%s' % deposition_id,
                params={'access_token': ZENODO_DEPOSIT}, data=json.dumps(data),
                headers=headers)
r.status_code
# 200


