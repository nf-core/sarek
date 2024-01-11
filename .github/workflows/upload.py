import requests
import os
import json

headers = {"Content-Type": "application/json"}
access_token = os.environ["ACCESS_TOKEN"]
params = {'access_token': access_token}
workspace_directory = os.environ["GITHUB_WORKSPACE"]
pipeline_version = os.environ["PIPELINE_VERSION"]
print("pipeline_version:" + pipeline_version)

# TODO: replace sandbox link https://zenodo.org/api/deposit/depositions
# https://sandbox.zenodo.org/api/deposit/depositions?access_token={access_token}
url = f"https://sandbox.zenodo.org/api/deposit/depositions"

# Create empty upload
r = requests.post(url,
                params=params,
                json={},
                headers=headers)

os.environ['DEPOSITION_ID'] = str(r.json()["id"])
deposition_id = r.json()["id"]

print("Create empty upload:\n")
print(r.json())
print("Deposition id: ")
print(deposition_id )

# Upload a new file
bucket_url = r.json()["links"]["bucket"]

filenames = [ #"deepvariant/HCC1395N/HCC1395N.deepvariant.vcf.gz",
                #"freebayes/HCC1395N/HCC1395N.freebayes.vcf.gz",
                #"haplotypecaller/HCC1395N/HCC1395N.haplotypecaller.filtered.vcf.gz",
                #"haplotypecaller/HCC1395N/HCC1395N.freebayes.vcf.gz",
                "strelka/HCC1395N/HCC1395N.strelka.genome.vcf.gz"]

for file in filenames:
    path = "./variant_calling/%s" % file
    with open(path, "rb") as fp:
        r = requests.put(
            "%s/%s" % (bucket_url, file),
            data=fp,
            params=params,
        )

# Add metadata to uploaded file
data = {
    'metadata': {
        'title': f'WES benchmark results nf-core/sarek v{pipeline_version}',
        'upload_type': 'data',
        'description': f'Variant calling results on benchmarking datasets produced with the nf-core/sarek v{pipeline_version}.',
        'creators': [{'name': 'Garcia, Maxime Ulysse',
                    'affiliation': 'Seqera, Barcelona'},
                    {'name': 'Hanssen, Friederike',
                    'affiliation': 'Quantitative Biology Center, Tuebingen'}]
    }
}

r = requests.put('https://sandbox.zenodo.org/api/deposit/depositions/%s' % deposition_id,
                params=params,
                data=json.dumps(data),
                headers=headers)

print("Add metadata: ")
print(r.status_code)
print(r.json())
print()

# filename = "HCC1395N.strelka.genome.vcf.gz"
# path = "./variant_calling/strelka/HCC1395N/%s" % filename
# print(os.listdir('../../'))
# print(os.listdir('../../../'))
# print(os.listdir('../'))
# print(os.listdir('./variant_calling/strelka/HCC1395N/'))

r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions/%s/actions/publish' % deposition_id,
                    params=params )

print("Publish data status code: ")
print(r.status_code)

