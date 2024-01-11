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
# Print DEPOSITION_ID to make it available for the workflow
print(f"::set-output name=deposition_id::{deposition_id}")

print("Create empty upload:\n")
print(r.json())
print("Deposition id: ")
print(deposition_id )

# Upload a new file
# bucket_url = r.json()["links"]["bucket"]

print(os.listdir('./variant_calling/strelka/HCC1395N/'))
# filenames = [ "deepvariant/NA12878_75M/NA12878_75M.deepvariant.vcf.gz",
#               "freebayes/NA12878_75M/NA12878_75M.freebayes.vcf.gz",
#               "haplotypecaller/NA12878_75M/NA12878_75M.haplotypecaller.filtered.vcf.gz",
#               "haplotypecaller/NA12878_75M/NA12878_75M.freebayes.vcf.gz",
#               "strelka/NA12878_75M/NA12878_75M.strelka.variants.vcf.gz",
#               "strelka/NA12878_75M/NA12878_75M.strelka.genome.vcf.gz",

#               "deepvariant/NA12878_200M/NA12878_200M.deepvariant.vcf.gz",
#               "freebayes/NA12878_200M/NA12878_200M.freebayes.vcf.gz",
#               "haplotypecaller/NA12878_200M/NA12878_200M.haplotypecaller.filtered.vcf.gz",
#               "haplotypecaller/NA12878_200M/NA12878_200M.freebayes.vcf.gz",
#               "strelka/NA12878_200M/NA12878_200M.strelka.variants.vcf.gz",
#               "strelka/NA12878_200M/NA12878_200M.strelka.genome.vcf.gz"]

# for file in filenames:
#     path = "./variant_calling/%s" % file
#     with open(path, "rb") as fp:
#         r = requests.put(
#             "%s/%s" % (bucket_url, file),
#             data=fp,
#             params=params,
#         )
#         print("%s/%s" % (bucket_url, file))
#         print(r.json())


# print("Upload new files")
# print(r.json())

# Add metadata to uploaded file

title = 'WES benchmark results nf-core/sarek v{}'.format(pipeline_version)
data = {
    'metadata': {
        'title': title,
        'upload_type': 'data',
        'description': 'Variant calling results on benchmarking datasets produced with nf-core/sarek',
        'creators': [{'name': 'Garcia, Maxime Ulysse', 'affiliation': 'Seqera, Barcelona'},
                    {'name': 'Hanssen, Friederike', 'affiliation': 'Quantitative Biology Center, Tuebingen'}]
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

# Publish this

# r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions/%s/actions/publish' % deposition_id,
#                     params=params )

# print("Publish data status code: ")
# print(r.status_code)

