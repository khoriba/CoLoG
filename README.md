![](https://github.com/khoriba/COLOG/blob/main/scripts/Colog.png)
# COLOG: Coronavirus Local Genomics analysis pipeline
CoLoG (Coronavirus Local Genomics) is a bioinformatics pipeline designed to perform basic SARS-CoV-2 genome analysis in local computing environments. It was developed to bring the genome-analysis component of Japan's COG-JP surveillance web service into an offline, site-operable workflow that laboratories can run securely on their own servers.

## Installation
On a Linux system with Conda installed, execute the following commands to install the CoLoG pipeline environment.

### 1. Clone the repository and build the Conda environment
```
$ git clone git@github.com:khoriba/COLOG.git
$ conda env create -n colog_env -f ./COLOG/environment.yml
```

### 2. Activate the environment and Run the setup scripts
```
$ conda activate colog_env
(colog_env)$ bash ./COLOG/setup.sh
```

### 3. Download the Nextclade reference dataset
If your system requires a proxy connection, include the --proxy option as shown below.
This will download the Nextclade reference dataset for the Wuhan-Hu-1 strain and place it under your Conda environment’s reference directory.
```
(colog_env)$ nextclade dataset get \
  --proxy 'http://proxy-server.co.jp:8080' \
  --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' \
  --output-dir "$CONDA_PREFIX/opt/reference/nextstrain/sars-cov-2/wuhan-hu-1/orfs"
```

## Usage
### 1. Preparing Input Data
Before running CoLoG, you need to organize your sequencing reads into individual directories by sample.
The following is an example command sequence that separates paired-end read data (gzip-compressed FASTQ files) and optional assembly data (FASTA files) into sample-specific directories. The `mkchank.sh` script automatically creates a subdirectory for each sample and moves the corresponding FASTQ files into that directory. After the script execution, each sample will have its own directory containing its paired-end files
```
# Check the input files
(colog_env)$ ls
JP01_R1_001.fastq.gz  JP02_R1_001.fastq.gz  JP03_R1_001.fastq.gz
JP01_R2_001.fastq.gz  JP02_R2_001.fastq.gz  JP03_R2_001.fastq.gz

# Run mkchank.sh
(colog_env)$ mkchank.sh

# Verify that the sample directories have been created
(colog_env)$ ls -R
.:
JP01  JP02  JP03

./JP01:
JP01_R1_001.fastq.gz  JP01_R2_001.fastq.gz

./JP02:
JP02_R1_001.fastq.gz  JP02_R2_001.fastq.gz

./JP03:
JP03_R1_001.fastq.gz  JP03_R2_001.fastq.gz
```

### 2. Example: Running the Pipeline with Paired-End Read Data
The following example demonstrates how to execute the CoLoG pipeline using paired-end read data (gzip-compressed FASTQ files).

In this example:
	•	The working directory contains a folder named RB01.
	•	The RB01 folder contains paired-end FASTQ files:
	•	Read 1 file: *_R1_*.fastq.gz
	•	Read 2 file: *_R2_*.fastq.gz
	•	The number of threads used is 16 (default is 4).
	•	The default primer set located at `$CONDA_PREFIX/opt/reference/nCovid-19-primerN7.mod.bed` is used for primer trimming.

```
# Check the input data
(colog_env)$ ls JP01
JP01_R1_001.fastq.gz  JP01_R2_001.fastq.gz

# Run the analysis pipeline
(colog_env)$ COLOG.sh --threads 16 JP01
```
Note:

If you have already organized your input data as described in Section 1, you can simply specify the folder containing each sample’s paired-end data to run the pipeline. CoLoG will automatically process the files within that directory.




