# Hybrid *de novo* assembly of transcriptome reference
## About

This repository contains a Nextflow pipeline for *de novo* transcriptome assembly based on long-read and short-read RNA-sequencing data. It includes various stages such as filtering, alignment, assembly, annotation, and polishing. The pipeline is implemented using Nextflow and Docker, and it is designed to be highly configurable and adaptable to different datasets and computational environments.

## System Requirements
-----
### Hardware Requirements

The exact memory requirement depends on the scale of dataset being processed. For the dataset we processed in the manuscript, the following specifications are recommended:

- RAM: 128+ GB
- CPU: 24+ threads

The runtime for the test dataset using a computer with 128GB RAM and 20 cores CPU (Intel 13600K) is approximately 3 mins.

### Software Requirements

Pipeline has been tested on the following operating system:
    - Linux: Fedora 36 (Linux kernel 5.17), 39 (Linux kernel 6.7.3)

- Nextflow (v23.04.1.5866+)
- Docker (v25.0.2+)

## Instructions for running this pipeline
-----

First clone this repository at a location of your choice:

```bash
git clone https://github.com/xuesoso/acoel_reference_assembly.git
cd acoel_reference_assembly/
```

To run the pipeline, you need to have [Docker](https://docs.docker.com/engine/install/) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) installed on your system. You can then build the Docker image using the provided script in the repository's root directory:

```bash
bash Docker/build_docker.sh
```

### Running the test pipeline

To test the pipeline, download the test dataset from this [google drive link](https://drive.google.com/file/d/11pL8G-3bgaXjzsJymLgYMJF-Xb8wx3X1/view?usp=sharing). Extract the **test_data.tar.gz** file and place the resulting folder under "data/" directory:

```bash
tar -xvf test_data.tar.gz
mv ${test_data} acoel_reference_assembly/data/
bash run_test.sh
```

### Running the pipeline on other dataset

To run the assembly and annotation workflow on other dataset, you need to generate a samplesheet for the fastq files (refer to `assets/test_samples.csv` as an example) and a yaml file that specifies fastq paths and conditions for Salmon grouper (i.e. `assets/grouper_config.yaml`). Please note that you may need to replace the paths for the input and output directories. Specify in `nextflow.config` the input directory for the dataset by modifying the `params.input_dir` variable.

#### Input samplesheet
The `assets/samples.csv` file specifies the paths for the input fastq files of Illumina, Nanopore, and Pacbio sequencing:

```csv
R1,R2,sample,platform
sample_1/R1.fastq.gz,sample_1/R2.fastq.gz,sample_1,illumina
sample_2/R1.fastq.gz,sample_2/R2.fastq.gz,sample_2,illumina
pacbio/RNA_seq.fastq.gz,,,pacbio
nanopore/RNA_seq.fastq.gz,,,nanopore
```

The `assets/grouper_config.yaml` file specifies the input paths of Illumina fastq for Salmon grouper, as well as the condition factors if there were different treatments for the samples. 

## Dataset
-----

Full dataset including the polished transcriptome reference, TransDecoder predicted peptide sequences, and the annotated functions for the products are deposited at NCBI GEO under accession number **[GSE242841](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242841)**.

## Citation
-----

A manuscript associated with this workflow is currently under review. If you find this repository useful, please consider citing our [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.06.21.545945v1):

**Coordinated wound responses in a regenerative animal-algal photosymbiotic metaorganism (2023).** Dania Nanes Sarfati, Yuan Xue, Eun Sun Song, Ashley Byrne, Daniel Le, Spyros Darmanis, Stephen R. Quake, Adrien Burlacot, James Sikes, Bo Wang. bioRxiv 2023.06.21.545945; doi: https://doi.org/10.1101/2023.06.21.545945

## License
-----

This project is licensed under the MIT license. For more details, please refer to the LICENSE file in the root directory of the repository.
