# Try the workflow

This small tutorial will explain to you how to run CAW on a small sample test data. Some variables are specific to Swedish UPPMAX cluster, but can be easily modified to suit any clusters.

## Make a test directory

```bash
mkdir test_CAW
cd test_CAW
```

## Test the workflow on a test tiny set

This workflow itself needs no installation. Nextflow will automatically fetch it from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

```bash
nextflow run SciLifeLab/CAW --test
```

If you're using a Swedish UPPMAX cluster, don't forget to provide your project ID.

```bash
nextflow run SciLifeLab/CAW --test --project <UPPMAX Project ID>
```

# Other possibilities for advanced users

## Clone the repository and test the workflow on a test tiny set

You can download the repository yourself from GitHub and run them directly:

```bash
git clone https://github.com/SciLifeLab/CAW test_CAW
cd test_CAW
nextflow run main.nf --test
```

## Using docker images

You can also use docker images:

```bash
nextflow run main.nf --test -profile docker
```

## Load Nextflow

If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:

```bash
module load Nextflow
```

Environnement variables are set up each time the module is loaded, so you might want to set them up after loading the module.

## Running tests in interactive mode on milou

You can try the test data by changing to the interactive mode on milou and run the test tiny set like:

```
$ interactive -A <UPPMAX Project ID> -p node
[ ... login messages ... ]
$ nextflow run main.nf -profile interactive --test --project <UPPMAX Project ID>
```

For more tests, see [Usage documentation](USAGE.md#test)
