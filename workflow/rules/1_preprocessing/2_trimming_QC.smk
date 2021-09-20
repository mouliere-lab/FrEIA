import pandas as pd

configfile: "../../../config/config.yaml" #Set config file.

ProjDirName = config["ProjName"] + "/trimmed"
tmp_dir = config["TmpDir"] + "/" + ProjDirName #Set TEMPDIR.
Samplesheet = pd.read_csv(config["Samplesheet"], delim_whitespace=True) #Read sample sheet in a dataframe.

localrules: all_trimming_quality_check

rule all_trimming_quality_check:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_report_interactive.html",
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_report_flat.html",
        expand(config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R{R_nr}_fastqc{ext}", zip,
            sample = Samplesheet["sample_name"],
            R_nr = [1,2]*len(Samplesheet["sample_name"]),
            ext = [".html",".zip"]*len(Samplesheet["sample_name"]))

rule fastqc_trimming:
    input:
        tmp_dir + "/1_trimming/{sample}_R1.fq.gz",
        tmp_dir + "/1_trimming/{sample}_R2.fq.gz"
    output:
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R1_fastqc.html",
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R1_fastqc.zip",
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R2_fastqc.html",
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R2_fastqc.zip"
    params:
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/"
    threads: config["ThreadNr"]
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}.log"

    shell:
        """
        fastqc {input} -t {threads} --outdir {params} \
        &> {log}
        """

rule multiqc_trimming:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R1_fastqc.html",
        sample=Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R1_fastqc.zip",
        sample=Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R2_fastqc.html",
        sample=Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/{sample}_R2_fastqc.zip",
        sample=Samplesheet["sample_name"])
    output:
        interactive = config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_report_interactive.html",
        flat = config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_report_flat.html"
    params:
        config["OutPath"] + "/" + ProjDirName + "/1_trimming_quality/fastqc/",
        config["OutPath"] + "/logs/" + ProjDirName + "/1_trimming_quality/fastqc/"
    resources:
        time=60
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc.tsv"
    log:
        interactive = config["OutPath"] + "/logs/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_interactive.log",
        flat = config["OutPath"] + "/logs/" + ProjDirName + "/1_trimming_quality/multiqc/multiqc_flat.log"
    shell:
        """
        multiqc {params} \
        --interactive \
        --force \
        -o $(dirname {output.interactive}) \
        --filename $(basename {output.interactive}) `# Next-level MultiQC hack` \
        --verbose &> {log.interactive} &

        multiqc {params} \
        --flat \
        --force \
        -o $(dirname {output.flat}) \
        --filename $(basename {output.flat}) `# Next-level MultiQC hack` \
        --verbose &> {log.flat} &
        wait
        """
