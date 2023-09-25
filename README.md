# Hybrid *de novo* assembly of transcriptome reference
## About

This repository contains a Nextflow pipeline for *de novo* transcriptome assembly based on long-read and short-read RNA-sequencing data. It includes various stages such as filtering, alignment, assembly, annotation, and polishing. The pipeline is implemented using Nextflow and Docker, and it is designed to be highly configurable and adaptable to different datasets and computational environments.

## Usage

First clone this repository at a location of your choice:

```bash
git clone https://github.com/xuesoso/acoel_reference_assembly.git
cd acoel_reference_assembly/
```

To test the pipeline, download the test dataset from this [google drive link](https://drive.google.com/file/d/11pL8G-3bgaXjzsJymLgYMJF-Xb8wx3X1/view?usp=sharing). Extract the **test_data.tar.gz** file and place the resulting folder in the "data" directory:

```bash
tar -xvf test_data.tar.gz
mv input data/
```

To run the pipeline, you need to have [Docker](https://docs.docker.com/engine/install/) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) installed on your system. You can then build the Docker image using the provided script in the repository's root directory:

```bash
bash Docker/build_docker.sh
```

After building the Docker image, you can run the test pipeline with the following command:

```bash
bash run_test.sh
```

To run the assembly and annotation workflow on other dataset, generate a spreadsheet of fastq locations (refer to `assets/test_samples.csv` as an example). Modify `assets/grouper_config.yaml` and add in the condition factors accordingly. Please note that you may need to replace the paths for the input and output directories.

## Dataset

Full dataset including the polished transcriptome reference, TransDecoder predicted peptide sequences, and the annotated functions for the products are deposited at NCBI GEO under accession number **[GSE242841](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242841)**.

## Citation

A manuscript associated with this workflow is currently under review. If you find this repository useful, please consider citing our [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.06.21.545945v1):

**Coordinated wound responses in a regenerative animal-algal photosymbiotic metaorganism (2023).** Dania Nanes Sarfati, Yuan Xue, Eun Sun Song, Ashley Byrne, Daniel Le, Spyros Darmanis, Stephen R. Quake, Adrien Burlacot, James Sikes, Bo Wang. bioRxiv 2023.06.21.545945; doi: https://doi.org/10.1101/2023.06.21.545945

## License

This project is licensed under the MIT license. For more details, please refer to the LICENSE file in the root directory of the repository.
