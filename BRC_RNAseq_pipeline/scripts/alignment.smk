import glob
import pathlib
import pandas as pd
import os
import time

os.environ["PATH"] = "/brcwork/bioinf/tools/pigz-2.4/:/brcwork/bioinf/tools/bcl2fastq-2.20/bin:" + os.environ["PATH"]

################################################################################
# helper functions

def extract_sample_info(samplesheet_path):
  samplesheet = pd.read_csv(samplesheet_path)
  data_line_num = samplesheet.index[samplesheet.iloc[:,0] == "[Data]"][0]
  sample_info = pd.read_csv(samplesheet_path,header = data_line_num + 2, dtype = "object")
  sample_info.fillna("", inplace=True)
  if "Lane" in list(sample_info):
    sample_info = sample_info.loc[[0]]
  return(sample_info)

def get_latest_run_name(archive_path, sample_info, sample):
  lab = sample_info.loc[sample,"Lab"]
  project = sample_info.loc[sample,"SampleProject"]
  fastq = "{}/Lab/{}/{}/fastq/{}/*".format(archive_path, lab, project, sample)
  fastq = glob.glob(fastq)
  assert len(fastq) > 0, "No fastq for sample {}".format(sample)
  fastq.sort()
  run_name = fastq[-1].split("/")[-1]
  return(run_name)

################################################################################
configfile: "config/config.yaml"

# setup

archive_path = config["archive"]

log_path = "{}/Log/{}".format(archive_path, config["log"])
os.makedirs(log_path, exist_ok = True)

samplesheet_path = config["samplesheet"]

# process samplesheet
sample_info = extract_sample_info(samplesheet_path)
samples = list(sample_info.loc[:,"Sample_ID"])
sample_info.index = samples

sample_info["Run"] = [get_latest_run_name(archive_path, sample_info, sample) for sample in samples]

################################################################################
# rules
localrules: all, complete

wildcard_constraints:
  sample = "[a-zA-Z0-9\-\_]+"

rule all:
  input:
    "{}/_complete".format(log_path)

rule merge_fastq:
  input:
    R1 = lambda wildcards: glob.glob("{}/Lab/{}/{}/fastq/{}/*/*_R1_001.fastq.gz".format(archive_path, wildcards.lab, wildcards.project, wildcards.sample)), 
    R2 = lambda wildcards: glob.glob("{}/Lab/{}/{}/fastq/{}/*/*_R2_001.fastq.gz".format(archive_path, wildcards.lab, wildcards.project, wildcards.sample))
  output:
    R1 = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/_fastq_tmp/R1.fastq.gz".format(archive_path),
    R2 = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/_fastq_tmp/R2.fastq.gz".format(archive_path)
  params:
    name = "mergefastq_{sample}",
    thread = "4",
    mem_per_thread = "1G",
    output = "{}/{{sample}}.merge_fastq.output".format(log_path),
    error = "{}/{{sample}}.merge_fastq.error".format(log_path)
  run:
    r1_fq = output.R1[:-3]
    shell("zcat {{input.R1}} > {}".format(r1_fq))
    shell("pigz -p4 {}".format(r1_fq))
    r2_fq = output.R2[:-3]
    shell("zcat {{input.R2}} > {}".format(r2_fq))
    shell("pigz -p4 {}".format(r2_fq))    

rule fastqc:
  input:
    R1 = rules.merge_fastq.output.R1,
    R2 = rules.merge_fastq.output.R2
  output:
    R1_report = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/fastqc/R1_fastqc.html".format(archive_path),
    R2_report = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/fastqc/R2_fastqc.html".format(archive_path)
  params:
    name = "fastqc_{sample}",
    thread = "1",
    mem_per_thread = "6G",
    output = "{}/{{sample}}.fastqc.output".format(log_path),
    error = "{}/{{sample}}.fastqc.error".format(log_path)
  run:
    out_dir = "/".join(output.R1_report.split("/")[:-1])
    shell("fastqc -o {} -f fastq --extract {{input}}".format(out_dir))

rule filter_fastq:
  input:
    R1 = rules.merge_fastq.output.R1,
    R2 = rules.merge_fastq.output.R2
  output:
    R1 = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/filtered/R1.fastq.gz".format(archive_path),
    R2 = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/filtered/R2.fastq.gz".format(archive_path)
  params:
    name = "filterfastq_{sample}",
    thread = "1",
    mem_per_thread = "6G",
    output = "{}/{{sample}}.filter_fastq.output".format(log_path),
    error = "{}/{{sample}}.filter_fastq.error".format(log_path)
  run:
    shell("bash scripts/filter.sh {input.R1} {input.R2} {output.R1} {output.R2}")

rule STAR:
  input:
    R1 = rules.filter_fastq.output.R1,
    R2 = rules.filter_fastq.output.R2
  output:
    bam = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/STAR_alignment/Aligned.out.bam".format(archive_path)
  params:
    name = "STAR_{sample}",
    thread = "8",
    mem_per_thread = "4G",
    output = "{}/{{sample}}.STAR.output".format(log_path),
    error = "{}/{{sample}}.STAR.error".format(log_path)
  run:
    ref_path = sample_info.loc[wildcards.sample,"Reference"] + "/star"
    out_prefix = "/".join(output.bam.split("/")[:-1]) + "/"
    shell("STAR --runThreadN 8 --genomeDir {} --readFilesIn {{input.R1}} {{input.R2}} --outFileNamePrefix {} --readFilesCommand gunzip -c --outSAMtype BAM Unsorted".format(ref_path, out_prefix))

rule htseq_count:
  input:
    bam = rules.STAR.output.bam
  output:
    count = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/htseq_count/counts.txt".format(archive_path)
  params:
    name = "htseqcount_{sample}",
    thread = "1",
    mem_per_thread = "4G",
    output = "{}/{{sample}}.htseq_count.output".format(log_path),
    error = "{}/{{sample}}.htseq_count.error".format(log_path)
  run:
    sample = wildcards.sample
    gtf_path = sample_info.loc[sample,"Reference"] + "/genes/genes.gtf"
    shell("bash scripts/htseqcount.sh {} {} {}".format(input.bam, gtf_path, output.count))

rule sort_bam:
  input:
    bam = rules.STAR.output.bam
  output:
    bam = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/STAR_alignment/possorted.bam".format(archive_path),
    bai = "{}/Lab/{{lab}}/{{project}}/results/{{sample}}/{{run_name}}/STAR_alignment/possorted.bam.bai".format(archive_path)
  params:
    name = "sortbam_{sample}",
    thread = "4",
    mem_per_thread = "4G",
    output = "{}/{{sample}}.sort_bam.output".format(log_path),
    error = "{}/{{sample}}.sort_bam.error".format(log_path)
  run:
    out_prefix = "/".join(output.bam.split("/")[:-1]) + "/_SORTtmp"
    os.makedirs(out_prefix, exist_ok = True)
    shell("samtools sort -@ {{params.thread}} -T {} -o {{output.bam}} -O bam {{input}}".format(out_prefix))
    shell("samtools index {output.bam}")

rule complete:
  input:
    fastqc = ["{}/Lab/{}/{}/results/{}/{}/fastqc/{}_fastqc.html".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, sample_info.loc[sample,"Run"], read) for sample in samples for read in ["R1", "R2"]], 
    counts = ["{}/Lab/{}/{}/results/{}/{}/htseq_count/counts.txt".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, sample_info.loc[sample,"Run"]) for sample in samples],
    bai = ["{}/Lab/{}/{}/results/{}/{}/STAR_alignment/possorted.bam.bai".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, sample_info.loc[sample,"Run"]) for sample in samples]
  output:
    "{}/_complete".format(log_path) 
  run:
    shell("echo $(date) > {output}")
