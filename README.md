# FrEIA - Fragment End Integrated Analysis

This workflow contains tools for pre-processing as well as
Fragment end sequence (FrES) analysis of cfDNA data.

1. Pre-processing of sequencing reads
1.1. Trimming by either BBduk or Cutadapt.
1.2. Mapping by bwa mem.
1.3. Reads are filtered by Samtools and Picard
1.4. Quality control by FastQC and Qualimap and is summarized by MultiQC.

2. FrES analysis
2.1. FrES analysis by the Fragment End Integrated Analysis (FrEIA) tool.

## Methodology and citation
Moldovan et al. Genome-wide cell-free DNA termini in patients with cancer. (2021)

## 1. Installing dependencies
Install [Anaconda][1].

## 2. Downloading the Core pipeline
First create a directory to store the workflow:
```mkdir FrEIA
cd FrEIA
Then download the workflow:
wget https://github.com/mouliere-lab/FrEIA.git
or
git clone https://github.com/mouliere-lab/FrEIA.git
```

## 3. Creating the conda environment for the workflow.
This can take one or two minutes...
`conda env create --name FrEIA --file /envs/environment.yaml`

See the [Anaconda website][2] for more help.

## 4. Activating the new environment:
`conda activate FrEIA`

## 5. Changing file extension to expected format
The FrEIA pipeline recognizes files with a `.fq.gz` extension. If your files
have a different extension use the script in `workflow/scripts/0_raw_seq/1_raw_renaming.sh`
to change it.

## 6. Setting the path to the files.
You can set the path to the files you want to use in `config/config.yaml`.
All the sample files must be in the same directory.

## 7. Creating a sample sheet
Create a new sample sheet.
The sample sheet is a space-delimited file with two labels in its header:
`group` `sample_name`

`group` represents the test group the given file belongs to. You can add as many
test groups as you prefer. Some analysis requires a control group.
Mark this with a `*` following the group name.

`sample_name` is the name of the unique name of the sample (file). For each pair of
paired end reads there is one unique sample name, which is the file name minus
the `_R1(2).fq.gz` suffix.

You can find an example in the `config` directory.
Set the path to the sample sheet in `config/config.yaml`.

## 8. Setting the reference genome you want to use.
You can set the path to the reference genome you want to use in `config/config.yaml`.

## 9. Setting the output location.
By default the results will be placed in the results directory.
You can change the path to the output directory in `config/config.yaml`.

## 10. Running preprocessing on a HPC cluster using Snakemake and slurm
Rules of the pipeline are run separately. To run a rule, first change directories
into the rule's folder
for example `cd workflow/rules/1_preprocessing/`
and run snakemake with the `-s` argument followed by the rule name
for example
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cluster "sbatch \
                     --time=60 \
                     --nodes=1 \
                     --partition=normal \
                     --mail-type=FAIL,TIME_LIMIT \
                     --mail-user=your@email.address \
                     --output=path/to/slurm_out/slurm-%j.out" \
          --jobs 500 \
          --max-jobs-per-second 3 \
          --max-status-checks-per-second 5 \
          -s workflow/rules/1_trimming.smk
```

## Inputs
The Core pipeline accepts paired `fq.gz` files.

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
        └── FrEIA/
            ├── 1_extract_fragment_ends/ - Data files containing data used by FrEIA.
            ├── 2_GC_bias/ - the GC bias of the sample groups.
            ├── 3_Abundances/
                ├── Genome_Lvl/ - genome level abundance data.
                ├── Chromosome_Lvl - chromosome level abundance data.
                └── SubChr_Lvl - subchromosomal abundance data.
            └── 4_Compare/ - the results of the FrEIA analysis comparing the groups of the data.
```
## Parameters and configuration
Parameters can be set in the `config.yaml` file found in the `config` directory.
To set or modify a path or parameter change the values after the `keyword:`.

## Bootstrap sampling
`compare.smk` offers a bootstrap sampling of the results of the preprocessing.
To do this, modify the `BsSampNr` parameter in the `config.yaml`.
`BsSampNr` default value is 0, meaning no bootstrap sampling. An integer above 0
results in multiple output files in equal numbers to the integer.
These outputs will be calculated from randomly sampled inputs, where the same read
can be represented multiple times (a.k.a. sampling with replacement).

## Calculate the FrEIA score
In the current version of the pipeline this step is done separately by
running the `5_FrEIA_score.py` script found at `workflow/scripts/FrEIA/`.
The script can be run using the following arguments:
`-ip` the path to the file containing the trinucleotide proportions found at
`results/trimmed/FrEIA/4_compare/`.
`id` the path to the file containing the diversity of the samples found at
`results/trimmed/FrEIA/4_compare/`.
`-st` the path to the sample table.
`-o` the path to the output file.
This version of the script uses a preset of trinucleotides for the score
calculation.

[1]: https://docs.anaconda.com/anaconda/install/linux/
[2]: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
