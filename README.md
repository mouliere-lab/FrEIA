# FrEIA - Fragment End Integrated Analysis

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE)
- [Citation](#citation)

## Overview
During cell death the DNA of tumor cells is released in to the bloodstream,
forming a pool of circulating tumor DNA (cfDNA). The cleavage of DNA originating
from cancer cells is aberrant compared to that of healthy cells, resulting in
shifted fragment end sequence signature. This can be leveraged to detect cancer
from low pass whole genome sequencing data.
The Fragment End Integrated Analysis tool extracts the proportion of fragment
end sequences and computes the diversity of these. Based on these metrics it
calculates a FrEIA score.

## Repo Contents
 - [FrEIA pipeline](./workflow)

## System Requirements

### Hardware Requirements
The FrEIA pipeline is designed to run on HPC clusters, as some steps can be
resource-heavy with real-world data sets.

If ran on a standard computer, we recommend:
RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

### Software Requirements

The package is tested on an Ubuntu 20.04.4 LTS operating system, but a it should
work on [Windows Subsystem for Linux (WSL)][1] as well.

## Installation Guide

### Package Dependencies
Prior to downloading FrEIA users should install [Anaconda][2].

FrEIA uses [Snakemake][3] for reproducibility.
Install Snakemake by running:
```
conda create -n snakemake -c bioconda snakemake
```
Activate your environment:
```
conda activate snakemake
```
Snakemake uses Mamba, a fast package manager. To install it type:
```
conda install mamba -n base -c conda-forge
```
Snakemake requires strict channel priorities in Conda, which is also good for
reproducibility. To set this, type:
```
conda config --set channel_priority strict
```
Every other dependency will be installed automatically during the first run of
the pipeline.
Software versions used by the pipeline can be found in the environment files
in the [envs](./workflow/envs/) folder.

In short the pipeline uses:
[BBduk][6] or [Cutadapt][7] for adapter trimming.
[BWA MEM][8] for mapping.
[Samtools][9] and [Sambamba][10] for quality filtering and marking duplicates.
[Qualimap][11] for quality control summarized by [MultiQC][12].

### Downloading FrEIA
Download FrEIA by:
```
wget https://github.com/mouliere-lab/FrEIA/archive/refs/heads/main.zip &&
unzip main.zip &&
mv FrEIA-main FrEIA &&
rm main.zip
```
or if you have Git installed:
```
git clone https://github.com/mouliere-lab/FrEIA.git
```

## Demo
The FrEIA pipeline accepts paired sequencing files with `.fq.gz` extension.
If your files have a different extension use the script in
`workflow/scripts/0_raw_seq/1_raw_renaming.sh` to change it.

### Dummy Data Set
To run this demo please download [these dummy files][4].
The dummy files are a synthetic admixtures of real sequencing files from healthy
controls and cancer patients respectively, downsampled to 100K reads.
The downloaded file contains a `metadata.csv` file used for the FrEIA score
calculation and machine learning classification.
This file contains the phenotype of the samples.

### Creating the Sample Sheet
The [sample sheet](./config/samplesheet.csv) is a space-delimited `.csv` file
with two labels in its header: `sample_name group`

An example containing the dummy file names can be found in the
[config](./config) directory.

`sample_name` is the unique name of the sample (file). For each pair of
input files there is one unique sample name, which is the file name minus
the `_R1(2).fq.gz` suffix.

`group` represents the single analysis group the given file belongs to. For the
final steps of the FrEIA analysis a control group is needed. Mark this with a
`*` following the group name.

### Using the Config File
The [config file](./config/config.yaml) is a `.yaml` format file containing the
parameters used during the run.

#### 1. Set the project name
This will be used as the name of the main folder containing results.
#### 2. Set the sample path
Add the path to the folder containing the `.fq.gz` files.
#### 3. Set the path to the sample sheet
Add the path to the sample sheet file.
#### 4. Set the path to the reference genome
Add the path to the indexed reference genome.
#### 5. Set the output path
Add the path to the output folder. A folder with the project name will be
created here containing all the results.
#### 6. Set the path to a temporary folder
The temporary folder (tmp) will be used to store temporary files. On HPC clusters
this should be on a `scratch drive` with fast read/write capability.
#### 7. Set the number of CPU cores the pipeline can use
The number of CPU cores should be less or equal to the available cores.
On a normal computer this is usually 4 or 8.
If you use a HPC cluster be aware that with the increase of the core number
the memory usage can also increase.
#### 8. Set the path to the metadata file
The metadata file is a space-delimited `.csv` file with two obligatory fields:
`sample_name phenotype`. The `sample_name` should be the same as in the
`samplesheet.csv` file. The `phenotype` can be anything, here we use control
and cancer. This file can contain batch information, for which the data will be
corrected for.
#### 9. Set the column containing the batch information
Fragment end sequence proportions are prone for batch effects. We correct this
using [pycombat][5]. You can set batch information in a new column in the
`metadata.csv` file. The name of this column should be set in the
[config file](./config/config.yaml).
Dummy files are coming from the same batch, so `None` is passed as parameter.

### Running FrEIA on a HPC cluster using Snakemake and Slurm
This is the recommended way of running the pipeline, as some steps are
resource-heavy.
Rules of the pipeline are run separately. To run a rule, first change directories
into the rule's folder for example `cd workflow/rules/1_preprocessing/`
and run snakemake with the `-s` argument followed by the rule name

For example trim the reads by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cluster "sbatch \
                     --time=60 \
                     --nodes=1 \
                     --partition=thin \
                     --cpus-per-task=32 \
                     --mem=64000 \
                     --mail-type=FAIL,TIME_LIMIT \
                     --mail-user=your@email.address \
                     --output=path/to/slurm_out/slurm-%j.out" \
          --jobs 500 \
          --max-jobs-per-second 3 \
          --max-status-checks-per-second 5 \
          -s 1_trimming.smk
```
### Running FrEIA on a local machine using Snakemake
Running on a local machine is not recommended as some steps may exceed available
resources.
Rules of the pipeline are run separately. To run a rule, first change
directories into the rule's folder and run snakemake with the `-s` argument
followed by the rule name.

In the Demo we are running FrEIA on a local machine.
Expected run time: ~ 1.5 h

#### 1. Trimming
Change directory: `cd workflow/rules/1_preprocessing/`

Remove sequencing adapters by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s 1_trimming.smk
```

#### 2. Mapping
Align the reads to the previously set reference genome by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s 2_mapping.smk
```

#### 3. Run FrEIA
Change directory: `cd ../2_FrEIA`.

Extract the fragment ends by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s 1_preprocessing.smk
```

Calculate fragment end proportions and the FrEIA score:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s 2_FrEIA.smk
```
### Results
In the output folder the following files and folders should be created:

```
.
└── output_folder/
    └── trimmed/
        ├── 1_mapping/ - Mapped bam files.
        ├── 1_mapping_flagstat/ - Flagstat output of mapped files.
        ├── 1_mapping_quality/ - QC results of the mapped files.
        │   └── multiqc/ - Aggregate of the QC results.
        └── 4_FrEIA/
            ├── 1_extract_fragment_ends/ - Data files containing fragment end data per read.
            ├── 3_Abundances/
                └── Genome_Lvl/ - Genome level aggregated fragment end sequence abundance data.
            ├── 4_Compare/
                └── Data/- Aggregated fragment end sequence proportions and diversity.
            └── 5_FrEIA_score/ - the FrEIA scores and related data (batch corrected trinucleotide end proportions, logFC of the trinucleotide end proportions).
```
Please download the expected results from [here][13]. The expected results
archive contains the expected folder structure, quality check results and the
final files with fragment end sequence proportions, diversity and FrEIA score
of the dummy data set.

## Citation
For usage of the pipeline or the FrEIA tool alone please cite:
Moldovan et al. Genome-wide cell-free DNA biological patterns in patients with cancer. (2022)

## Machine learning classification
Classification is performed separately from the workflow by running
`2_classification/classify.py`.

Before running the script install dependencies:
`conda env create -f /envs/Classify_env.yaml -n Classify`
`conda activate Classify`

To run the classifier:
`python3 workflow/scripts/2_classification/classify.py \
    -i <path/to/input.csv> \
    -o <path/to/output_folder/file_name_prefix_> \
    -m RF \
    -y <Column name of the output variables>
    -ign <Name of the column to ignore, usually sample_name>`

The input file should contain at least 3 columns:
`sample_name`: the sample IDs.
`output_vars`: the output variables (y). These should be values of 0 or 1.
(ex. control and cancer).
`input_vars`: the input variables. These can be multiple columns.
(ex. FrEIA_score, ichorCNA_TF, P20_150, etc.)

The selected estimator above is the Random Forest estimator `-m RF`.
The estimator and parameters were chosen based on performance of the models as
described in the Methods section of Moldovan et al. 2022.

The outputs of the classification:
`_confmat.csv`: the confusion matrix of the prediction.
`_metrics_predict.csv`: different performance metrics of the prediction.
`_roc_predict.csv`: the FPR and TPR rates at different thresholds.

[1]: https://docs.microsoft.com/en-us/windows/wsl/install
[2]: https://docs.anaconda.com/anaconda/install/linux/
[3]: https://snakemake.github.io/
[4]: https://doi.org/10.6084/m9.figshare.20090063.v1
[5]: https://github.com/epigenelabs/pyComBat
[6]: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/
[7]: https://cutadapt.readthedocs.io/en/stable/
[8]: http://bio-bwa.sourceforge.net/bwa.shtml
[9]: http://www.htslib.org/
[10]: https://lomereiter.github.io/sambamba/
[11]: http://qualimap.conesalab.org/
[12]: https://multiqc.info/
[13]: https://doi.org/10.6084/m9.figshare.20090162.v1
