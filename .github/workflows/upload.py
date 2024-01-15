import requests
import os
import json

headers = {"Content-Type": "application/json"}
access_token = os.environ["ACCESS_TOKEN"]
params = {"access_token": access_token}
workspace_directory = os.environ["GITHUB_WORKSPACE"]
pipeline_version = os.environ["PIPELINE_VERSION"]

# TODO: replace sandbox link f"https://zenodo.org/api/deposit/depositions"
url = f"https://sandbox.zenodo.org/api/deposit/depositions"

# Create empty upload
r = requests.post(url, params=params, json={}, headers=headers)

print("Create empty upload:\n")
print(r.json())
print(r.status_code)

deposition_id = r.json()["id"]

## Store deposition ID
with open("deposition_id.txt", "w") as f:
    f.write(str(deposition_id))

# Upload a new file
bucket_url = r.json()["links"]["bucket"]

filenames = [
    "deepvariant/NA12878_75M/NA12878_75M.deepvariant.vcf.gz",
    "freebayes/NA12878_75M/NA12878_75M.freebayes.vcf.gz",
    "haplotypecaller/NA12878_75M/NA12878_75M.haplotypecaller.filtered.vcf.gz",
    "strelka/NA12878_75M/NA12878_75M.strelka.variants.vcf.gz",
    "deepvariant/NA12878_200M/NA12878_200M.deepvariant.vcf.gz",
    "freebayes/NA12878_200M/NA12878_200M.freebayes.vcf.gz",
    "haplotypecaller/NA12878_200M/NA12878_200M.haplotypecaller.filtered.vcf.gz",
    "strelka/NA12878_200M/NA12878_200M.strelka.variants.vcf.gz",
]

for file in filenames:
    path = "./variant_calling/%s" % file
    with open(path, "rb") as fp:
        r = requests.put(
            "%s/%s" % (bucket_url, os.path.basename(file)),
            data=fp,
            params=params,
        )
        print(r.json())

# Add metadata to uploaded file
title = "WES benchmark results nf-core/sarek v{}".format(pipeline_version)
data = {
    "metadata": {
        "title": title,
        "upload_type": "dataset",
        "description": "Variant calling results on benchmarking datasets produced with nf-core/sarek",
        "creators": [
            {"name": "Garcia, Maxime Ulysse", "affiliation": "Seqera, Barcelona"},
            {"name": "Hanssen, Friederike", "affiliation": "Quantitative Biology Center, Tuebingen"},
        ],
    }
}

# TODO replace sandbox link https://zenodo.org/api/deposit/depositions/
r = requests.put(
    "https://sandbox.zenodo.org/api/deposit/depositions/%s" % deposition_id,
    params=params,
    data=json.dumps(data),
    headers=headers,
)

print("Add metadata: ")
print(r.status_code)
print(r.json())

# TODO only uncomment once everything works, replace sandbox link
# Publish this
r = requests.post(
    "https://sandbox.zenodo.org/api/deposit/depositions/%s/actions/publish" % deposition_id, params=params
)

print("Publish data status code: ")
print(r.status_code)
