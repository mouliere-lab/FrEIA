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
        #-m $((({resources.mem_mb}-8000)/{threads}))M \
rule mark_duplicates:
    input:
        tmp_dir + "/1_mapping/{sample}.bam"
    output:
        bam = config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
        mtrx = config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}_dup_metrics.txt"
    priority: 100
    threads: config["ThreadNr"]
    resources: max_file_handles=50000, # SURFsara interactive: ulimit -n = 51200
               mem_mb=54000
    params: jar = "$CONDA_PREFIX/share/picard-2.22.2-0/picard.jar",
            TMP_DIR = tmp_dir + "/2_mark_duplicates/"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/2_mark_duplicates/{sample}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log"
    shell:
        """
        java -Xmx{resources.mem_mb}M\
        -XX:ParallelGCThreads={threads}\
        -jar {params.jar} MarkDuplicates\
        INPUT={input}\
        REMOVE_DUPLICATES=false\
        MAX_FILE_HANDLES={resources.max_file_handles}\
        CREATE_INDEX=false\
        COMPRESSION_LEVEL=9\
        TMP_DIR={params.TMP_DIR}\
        METRICS_FILE={output.mtrx}\
        OUTPUT={output.bam}\
        &> {log}
        """

rule flagstat:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping_flagstat/{sample}_flagstats.txt"
    threads: config["ThreadNr"]
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
        sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log",
        sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}_dup_metrics.txt",
        sample = Samplesheet["sample_name"])
    output:
        interactive = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_interactive.html",
        flat = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_flat.html"
    resources:
        time=60
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
