import requests
import os

headers = {"Content-Type": "application/json"}
params = {'access_token': os.environ["ACCESS_TOKEN"]}

r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions',
                params=params,
                json={},
                headers=headers)
r.status_code
# 201
bucket_url = r.json()["links"]["bucket"]

filename = "*.vcf.gz"
path = "./variant_calling/%s" % filename

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
                params={'access_token': ACCESS_TOKEN}, data=json.dumps(data),
                headers=headers)
r.status_code
# 200


