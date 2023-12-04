import requests
import os

headers = {"Content-Type": "application/json"}
access_token = os.environ["ACCESS_TOKEN"]
workspace_directory = os.environ["GITHUB_WORKSPACE"]

# TODO: replace sandbox link https://zenodo.org/api/deposit/depositions
# https://sandbox.zenodo.org/api/deposit/depositions?access_token={access_token}
url = f"https://sandbox.zenodo.org/api/deposit/depositions/"

# Create empty upload
params = {'access_token': access_token}

r = requests.post(url,
                    params=params,
                    json={},
                    headers=headers)

# Upload a new file
bucket_url = r.json()["links"]["bucket"]

filename = "HCC1395N.strelka.genome.vcf.gz"
path = "./variant_calling/strelka/HCC1395N/%s" % filename
with open(path, "rb") as fp:
    r = requests.put(
        "%s/%s" % (bucket_url, filename),
        data=fp,
        params=params,
    )
r.json()

# Add metadata to uploaded file

data = {
    'metadata': {
        'title': 'My first upload',
        'upload_type': 'poster',
        'description': 'This is my first upload',
        'creators': [{'name': 'Doe, John',
                    'affiliation': 'Zenodo'}]
    }
}

r = requests.put('https://sandbox.zenodo.org/api/deposit/depositions/%s' % deposition_id,
                params=params, data=json.dumps(data),
                headers=headers)
r.status_code


# filename = "HCC1395N.strelka.genome.vcf.gz"
# path = "./variant_calling/strelka/HCC1395N/%s" % filename
# print(os.listdir('../../'))
# print(os.listdir('../../../'))
# print(os.listdir('../'))
# print(os.listdir('./variant_calling/strelka/HCC1395N/'))

os.environ['DEPOSITION_ID'] = r.json()["metadata"]["recid"]
# TODO add publication step
#r = requests.post('https://zenodo.org/api/deposit/depositions/%s/actions/publish' % deposition_id,
#                      params={'access_token': ACCESS_TOKEN} )
#r.status_code
# 202

