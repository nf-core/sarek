---
order: 1
---

# Variant Calling Tutorial

These pages are a tutorial workshop for the [Nextflow](https://www.nextflow.io) pipeline [nf-core/sarek](https://nf-co.re/sarek).

In this workshop, we will recap the application of next generation sequencing to identify genetic variations in a genome. You will learn how to use the pipeline sarek to carry out this data-intensive workflow efficiently. We will cover topics such as experimental design, configuration of the pipeline and code execution.

Please note that this is not an introductory workshop, and we will assume some basic familiarity with Nextflow.

By the end of this workshop, you will be able to:

- understand the key concepts behind variant calling, as adopted in this pipeline
- analyse simple NGS datasets with the sarek workflow
- customise some of its features for your own variant calling analyses
- integrate different sources of information to identify candidate variants
- make a hypothesis about variant interpretation using the output of sarek

Let's get started!

## Running with Gitpod

In order to run this using GitPod, please make sure:

1. You have a GitHub account: if not, create one [here](https://github.com/signup)
2. Once you have a GitHub account, sign up for GitPod using your GitHub user [here](https://gitpod.io/login/) choosing "continue with GitHub".

Now you're all set and can use the following button to launch the service:

[![Open in GitPod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/lescai-teaching/sarek-tutorial)

## Additional documentation

- You can find detailed documentation on **Nextflow** [here](https://www.nextflow.io/docs/latest/)
- You can find additional training on [these pages](https://training.nextflow.io)

## Credits & Copyright

This training material has been written by [Francesco Lescai](https://github.com/lescai) during the [nf-core](https://nf-co.re) Hackathon in Barcelona, 2023. It was originally meant as a contribution for the nf-core community, and aimed at anyone who is interested in using nf-core pipelines for their studies or research activities.

The Docker image and Gitpod environment used in this repository have been created by [Seqera](https://seqera.io) but have been made open-source ([CC BY-NC-ND](https://creativecommons.org/licenses/by-nc-nd/4.0/)) for the community.

All examples and descriptions are licensed under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).
