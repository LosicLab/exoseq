# LosicLab/exoseq

## Introduction

**LosicLab/ExoSeq** is a bioinformatics analysis pipeline that performs best-practice analysis pipeline for Exome Sequencing data. It is forked from nfcore/ExoSeq.

The pipeline is built based on [GATK](https://software.broadinstitute.org/gatk/best-practices/) best practices using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. The main steps done by pipeline are the following (more information about the processes can be found [here](docs/processes.md)).

* Alignment - bwa
* Marking Duplicates - picard
* Recalibration - gatk 4
* Realignment - gatk 4
* Variant Calling (Somatic or SNP) - gatk 4
* Variant Filtration - gatk 4

## Documentation
The LosicLab pipeline comes with the documentation forked from the original nf-core repository, found in the `docs/` directory:

1. [Pipeline installation and configuration instructions](docs/installation.md)
2. Pipeline configuration
   * [Local installation](docs/configuration/local.md)
   * [Amazon Web Services](docs/configuration/aws.md)
   * [Swedish UPPMAX clusters](docs/configuration/uppmax.md)
   * [Swedish cs3e Hebbe cluster](docs/configuration/c3se.md)
   * [Tübingen QBiC clusters](docs/configuration/qbic.md)
   * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
   * [Preparing custom exome capture kits](docs/kits.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

The pipeline now also has support for the MSSM Minerva HPC.

## Credits
The original nf-core/exoseq pipeline was initally developed by Senthilkumar Panneerselvam ([@senthil10](https://github.com/senthil10)) with a little help from Phil Ewels ([@ewels](https://github.com/ewels)) at the National Genomics Infrastructure, part of SciLifeLab in Stockholm and has been extended by Alex Peltzer ([@apeltzer](https://github.com/apeltzer)), Marie Gauder ([@mgauder](https://github.com/mgauder)) from QBIC Tuebingen/Germany as well as Marc Hoeppner ([@marchoeppner](https://github.com/marchoeppner)) from IKMB Kiel/Germany.

Many thanks also to others who have helped out along the way too, including [@pditommaso](https://github.com/pditommaso), [@colindaven](https://github.com/colindaven).
