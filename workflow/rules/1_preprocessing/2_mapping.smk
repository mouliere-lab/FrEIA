import pandas as pd

configfile: "../../../config/config.yaml" #Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], delim_whitespace=True) #Read sample sheet in a dataframe.

localrules: all_mapping

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed"
    tmp_dir = config["TmpDir"] + "/" + ProjDirName #Set TEMPDIR.
    MapIn = tmp_dir + "/1_trimming/"

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed"
    tmp_dir = config["TmpDir"] + "/" + ProjDirName #Set TEMPDIR.
    MapIn = config["SamplePath"] + "/"

else:
    print("\nMissing trimming option in config file!\n")
    exit()

rule all_mapping:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_interactive.html",
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
                sample = Samplesheet["sample_name"])

rule map:
    input:
        MapIn + "{sample}_R1.fq.gz",
        MapIn + "{sample}_R2.fq.gz"
    output:
        tmp_dir + "/1_mapping/{sample}.bam"
    params:
        ref=config["RefPath"]
    threads: config["ThreadNr"]
    conda: "../../envs/Preprocessing_env.yaml"
    resources: mem_mb=36000+8000
               #> Prevent memory exhaustion by samtools sort
               #>> 1500MB/thread for sorting, 8000 for BWA
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping/{sample}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping/{sample}.log"
    shell:
        """
        bwa mem {params.ref} {input} \
        -t {threads} \
        2> {log} | \
        samtools sort \
        -@ {threads} \
        -o {output} 2>> {log}
        """


rule mark_duplicates:
    input:
        tmp_dir + "/1_mapping/{sample}.bam"
    output:
        bam = config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
    threads: config["ThreadNr"]
    params:
        tmpDir = tmp_dir + "/2_mark_duplicates"
    conda: "../../envs/Preprocessing_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/2_mark_duplicates/{sample}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log"
    shell:
        """
        sambamba markdup \
        -t {threads} \
        --tmpdir {params.tmpDir} \
        {input} \
        {output} \
        2>> {log}
        """

rule flagstat:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping_flagstat/{sample}_flagstats.txt"
    threads: config["ThreadNr"]
    conda: "../../envs/Preprocessing_env.yaml"
    resources:
        time = 10
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_flagstat/{sample}.tsv"
    shell:
        """
        samtools flagstat {input} --threads {threads} > {output}
        """

rule qualimap:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        dir=directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}")
    threads: config["ThreadNr"]
    conda: "../../envs/Preprocessing_env.yaml"
    resources:
            mem_mb=60000, #round(60/(math.sqrt(12))),
            time=120
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}.log"
    shell:
        """
        unset DISPLAY; \
        qualimap bamqc -bam {input} \
        --java-mem-size={resources.mem_mb}M \
        -nt {threads} \
        -outdir {output.dir} \
        `# -outfile {output} #> gives errors` \
        &> {log}
        """

rule multiqc_mapping:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}",
        sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_flagstat/{sample}_flagstats.txt",
        sample = Samplesheet["sample_name"])
    output:
        interactive = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_interactive.html",
        flat = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_flat.html"
    resources:
        time=60
    conda: "../../envs/Preprocessing_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_multiqc/multiqc.tsv"
    log:
        interactive = config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_multiqc/multiqc_interactive.log",
        flat = config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_multiqc/multiqc_flat.log"
    shell:
        """
        multiqc {input} \
        --interactive \
        --force \
        -o $(dirname {output.interactive}) \
        --filename $(basename {output.interactive}) `# Next-level MultiQC hack` \
        --verbose &> {log.interactive} &

        multiqc {input} \
        --flat \
        --force \
        -o $(dirname {output.flat}) \
        --filename $(basename {output.flat}) `# Next-level MultiQC hack` \
        --verbose &> {log.flat} &
        wait
        """
