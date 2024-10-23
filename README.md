[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**emilytrybulec/repeat_curation** is a bioinformatics pipeline that takes a finished genome and performs repeat analysis. It produces a masked genome (.fasta), files containing coordinates of regions identified as repeats (.bed) for further manual curation, and images depicting output from TE Trimmer (.pdf). 

1. [`Repeat Modeler BuildDatabase`](https://github.com/Dfam-consortium/RepeatModeler/tree/master?tab=readme-ov-file#example-run) 
2. [`Repeat Modeler`](https://github.com/Dfam-consortium/RepeatModeler)
3. [`TE Trimmer`](https://github.com/qjiangzhao/TEtrimmer)
4. [`Repeat Masker`](https://www.repeatmasker.org/)
   

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

First, go through nextflow.config to configure the pipeline to your needs. Each option can be modified to change which programs run and command line options.  


`nextflow.config`:

```config
params {

    // Input options
    te_trimmer                 = false
    repeat_masker              = true
    cons_thr                   = 0.5

    soft_mask                  = true

    species                    = null
    genome_fasta               = null
    consensus_fasta            = null
}
```


Next, create a params.yaml file to input information in place of the null configurations. This file will, at a minimum, contain your genome and preferred out directory name. Optionally, a consensus path and RepeatMasker species flag can be supplied.    
  
`params.yaml`:

```yaml
params {
   genome_fasta           : "/core/projects/colossalanalyses/Finished_Genomes_for_Annotation/BayDuikerCDO11_5Jan2023_RaconR3.fasta"
   outdir                 :  "bay_duiker_softmask"
   species                : "cow"
}
```

Now, you can run the pipeline using:

```bash
nextflow pull emilytrybulec/repeat_curation
nextflow run emilytrybulec/repeat_curation \
   -profile <docker/singularity/.../institute> \
   -params-file params.yaml
```

Xanadu users: please refer to the [`example script`](https://github.com/emilytrybulec/repeat_curation/blob/main/nextflow.sh).

## Credits

emilytrybulec/repeat_curation was originally written by Emily Trybulec with the help of Jessica Storer.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
If you use emilytrybulec/repeat_curation for your analysis, please cite it using this git.


This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
