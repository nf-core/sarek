# Variant calling <!-- omit in toc -->

- [Germline variant calling](#germline-variant-calling)
- [Somatic variant calling with tumor - normal pairs](#somatic-variant-calling-with-tumor---normal-pairs)
- [Somatic variant calling with tumor only samples](#somatic-variant-calling-with-tumor-only-samples)

## Germline variant calling

Using Sarek, germline variant calling will always be performed if a variant calling tool with a germline mode is selected.
You can specify the variant caller to use with the `--tools` parameter (see [usage](./usage.md) for more information).

Germline variant calling can currently only be performed with the following variant callers:

- FreeBayes
- HaplotypeCaller
- Manta
- mpileup
- Sentieon (check the specific [sentieon](sentieon.md) documentation)
- Strelka
- TIDDIT

For more information on the individual variant callers, and where to find the variant calling results, check the [output](output.md) documentation.

## Somatic variant calling with tumor - normal pairs

Using Sarek, somatic variant calling will be performed, if your input tsv file contains tumor / normal pairs (see [input](input.md) documentation for more information).
Different samples belonging to the same patient, where at least one is marked as normal (`0` in the `Status` column) and at least one is marked as tumor (`1` in the `Status` column) are treated as tumor / normal pairs.

If tumor-normal pairs are provided, both germline variant calling and somatic variant calling will be performed, provided that the selected variant caller allows for it.
If the selected variant caller allows only for somatic variant calling, then only somatic variant calling results will be generated.

Here is a list of the variant calling tools that support somatic variant calling:

- ASCAT (check the specific [ASCAT](ascat.md) documentation)
- ControlFREEC
- FreeBayes
- Manta
- MSIsensor
- Mutect2
- Platypus
- Sentieon (check the specific [sentieon](sentieon.md) documentation)
- Strelka

For more information on the individual variant callers, and where to find the variant calling results, check the [output](output.md) documentation.

## Somatic variant calling with tumor only samples

Somatic variant calling with only tumor samples (no matching normal sample), is not recommended by the GATK best practices.
This is just supported for a limited variant callers.

Here is a list of the variant calling tools that support tumor-only somatic variant calling:

- Manta
- mpileup
- Mutect2
- Platpus
- TIDDIT

For more information on the individual variant callers, and where to find the variant calling results, check the [output](output.md) documentation.
