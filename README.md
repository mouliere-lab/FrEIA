# FrEIA - Fragment End Integrated Analysis

This workflow contains tools for pre-processing as well as
Fragment end sequence analysis of cfDNA data.

1. Pre-processing of sequencing reads
1.1. Trimming by either BBduk or Cutadapt.
1.2. Mapping by bwa mem.
1.3. Reads are filtered by Samtools and Picard
1.4. Quality control by FastQC and Qualimap and is summarized by MultiQC.

2. Fragment end sequence analysis
2.1. Fragment end sequence analysis by the Fragment End Integrated Analysis (FrEIA) tool.
2.2. FrEIA score calculation

3. Classification toolkit

## Methodology and citation
Moldovan et al. Genome-wide cell-free DNA termini in patients with cancer. (2021)

## 1. Installing dependencies
Install [Anaconda][1].

## 2. Downloading the FrEIA pipeline
First create a directory to store the workflow:
```mkdir FrEIA
cd FrEIA
Then download the workflow:
wget https://github.com/mouliere-lab/FrEIA.git
or
git clone https://github.com/mouliere-lab/FrEIA.git
```
## 3. Changing input file extensions to expected format
The FrEIA pipeline recognizes files with a `.fq.gz` extension. If your files
have a different extension use the script in `workflow/scripts/0_raw_seq/1_raw_renaming.sh`
to change it.

## 4. Setting the path to the input files
You can set the path to the files you want to use in `config/config.yaml`.
All the input files must be in the same directory.

## 5. Creating a sample sheet
The sample sheet is a space-delimited file with two labels in its header:
`group` `sample_name`

An example can be found in the `config` directory.

`group` represents the single test group the given file belongs to. You can add as many
test groups as you prefer. For the final step of the FrEIA analysis a control
group is needed. Mark this with a `*` following the group name.

`sample_name` is the name of the unique name of the sample (file). For each pair of
input files there is one unique sample name, which is the file name minus
the `_R1(2).fq.gz` suffix.

Set the path to the sample sheet in `config/config.yaml`.

## 6. Setting the reference genome
You can set the path to the reference genome in `config/config.yaml`.

## 7. Setting the output location.
By default the results will be placed in the results directory.
You can change the path to the output directory in `config/config.yaml`.

## 8. Running the FrEIA pipeline on a HPC cluster using Snakemake and Slurm
This is the recomended way of running the pipeline, as some steps are resource-heavy.
Rules of the pipeline are run separately. To run a rule, first change directories
into the rule's folder
for example `cd workflow/rules/1_preprocessing/`
and run snakemake with the `-s` argument followed by the rule name
for example trim the reads by running:
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
## 9. Running the FrEIA pipeline on a local machine using Snakemake and Slurm
Running on a local machine is not recomended as some steps may exceed ther
available resources.
Rules of the pipeline are run separately. To run a rule, first change directories
into the rule's folder
for example `cd workflow/rules/1_preprocessing/`
and run snakemake with the `-s` argument followed by the rule name
for example trim the reads by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s 1_trimming.smk
```

## Inputs
The FrEIA pipeline accepts paired `fq.gz` files.

## Outputs
The workflow creates the following output structure and outputs:
```
.
└── results/
    └── trimmed/
        ├── 1_trimming_quality/ - QC results of the trimmed files.
        │   └── multiqc/ - aggregate of the QC results.
        ├── 1_mapping/ - Mapped bam files.
        ├── 1_mapping_flagstat/ - Flagstat output of mapped files.
        ├── 1_mapping_quality/ - QC results of the mapped files.
        │   └── multiqc/ - aggregate of the QC results.
        └── 4_FrEIA/
            ├── 1_extract_fragment_ends/ - Data files containing data used by FrEIA.
            ├── 3_Abundances/
                └── Genome_Lvl/ - genome level abundance data.
            ├── 4_Compare/ - the results of the FrEIA analysis comparing the groups of the data.
            └── 5_FrEIA_score/ - the FrEIA scores and related data (batch corrected trinucleotide end proportions, logFC of the trinucleotide end proportions).
```
## Parameters and configuration
Parameters can be set in the `config.yaml` file found in the `config` directory.
To set or modify a path or parameter change the values after the `keyword:`.

# Machine learning classification
Classification can be performed separately from the workflow by running
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

[1]: https://docs.anaconda.com/anaconda/install/linux/
[2]: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
