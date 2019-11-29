import glob
import pathlib
import pandas as pd
import os

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

################################################################################
# setup

run_path = os.path.abspath(config["run"])
run_name = run_path.split("/")[-1]

# paths
archive_path = "/brcwork/sequence/Archive"

bcl2fastq_out_path = "{}/Run/{}/bcl2fastq".format(archive_path, run_name)
samplesheet_path = config["samplesheet"]

# process samplesheet
sample_info = extract_sample_info(samplesheet_path)
samples = list(sample_info.loc[:,"Sample_ID"])
sample_info.index = samples

################################################################################
# rules
wildcard_constraints:
  sample = "[a-zA-Z0-9\-\_]+"

rule all:
  input:
    counts = ["{}/Lab/{}/{}/results/{}/{}/htseq_count/counts.txt".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, run_name) for sample in samples],
    fastqc = ["{}/Lab/{}/{}/results/{}/{}/fastqc/summary.txt".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, run_name) for sample in samples],
    bai = ["{}/Lab/{}/{}/results/{}/{}/STAR_alignment/possorted.bam.bai".format(archive_path, sample_info.loc[sample,"Lab"], sample_info.loc[sample,"SampleProject"], sample, run_name) for sample in samples]


rule bcl2fastq:
  input:
    samplesheet_path,
    run_path + "/RTAComplete.txt"
  output:
    bcl2fastq_out_path + "/Stats/Stats.json"
  params:
    thread = "6",
    mem_per_thread = "6G"
  run:
    print({input})
    
    shell("bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --auto-set-to-zero-barcode-mismatches --find-adapters-with-sliding-window --adapter-stringency 0.9 --mask-short-adapter-reads 35 --minimum-trimmed-read-length 35 -R {} --sample-sheet {} -o {}".format(run_path,samplesheet_path,bcl2fastq_out_path))
    
    for sample in samples:
      print(sample)
      project = sample_info.loc[sample,"SampleProject"]
      lab = sample_info.loc[sample,"Lab"]
      for fq in glob.glob("{}/{}_S*_001.fastq.gz".format(bcl2fastq_out_path, sample)):
        print(fq)
        fq_file = fq.split("/")[-1]
        dest = "{}/Lab/{}/{}/fastq/{}/{}".format(archive_path, lab, project, sample, run_name)
        os.makedirs(dest, exist_ok = True)
        os.rename(fq, "{}/{}".format(dest,fq_file))

rule fastqc:
  input:
    rules.bcl2fastq.output
  output:
    archive_path + "/Lab/{Lab}/{project}/results/{sample}/" + run_name + "/fastqc/summary.txt"
  params:
    thread = "1",
    mem_per_thread = "6G"
  run:
    shell("touch {output}")

rule STAR:
  input:
    rules.bcl2fastq.output
  output:
    bam = archive_path + "/Lab/{Lab}/{project}/results/{sample}/" + run_name + "/STAR_alignment/Aligned.out.bam"
  params:
    thread = "8",
    mem_per_thread = "4G"
  run:
    print(input)
    
    lab = wildcards.Lab
    project = wildcards.project
    sample = wildcards.sample
    ref_path = sample_info.loc[sample,"Reference"] + "/star"
    fastq_path = "{}/Lab/{}/{}/fastq/{}".format(archive_path, lab, project, sample)
    out_prefix = "/".join(output.bam.split("/")[:-1]) + "/" 
    shell("bash scripts/STAR.sh {} {} {}".format(fastq_path, out_prefix, ref_path))

rule htseq_count:
  input:
    bam = rules.STAR.output.bam
  output:
    count = archive_path + "/Lab/{Lab}/{project}/results/{sample}/" + run_name + "/htseq_count/counts.txt"
  params:
    thread = "1",
    mem_per_thread = "4G"
  run:
    sample = wildcards.sample
    gtf_path = sample_info.loc[sample,"Reference"] + "/genes/genes.gtf"
    shell("bash scripts/htseqcount.sh {} {} {}".format(input.bam, gtf_path, output.count))

rule sort_bam:
  input:
    bam = rules.STAR.output.bam
  output:
    bam = archive_path + "/Lab/{Lab}/{project}/results/{sample}/" + run_name + "/STAR_alignment/possorted.bam"
  params:
    thread = "4",
    mem_per_thread = "4G"
  run:
    out_prefix = "/".join(output.bam.split("/")[:-1]) + "/_SORTtmp"
    os.makedirs(out_prefix, exist_ok = True)
    shell("samtools sort -@ {{params.thread}} -T {} -o {{output}} -O bam {{input}}".format(out_prefix))

rule index_bam:
  input:
    bam = rules.sort_bam.output.bam
  output:
    bai = archive_path + "/Lab/{Lab}/{project}/results/{sample}/" + run_name + "/STAR_alignment/possorted.bam.bai"
  params:
    thread = "1",
    mem_per_thread = "2G"
  run:
    shell("samtools index {input}")    


