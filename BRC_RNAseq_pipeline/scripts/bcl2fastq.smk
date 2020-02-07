import glob
import pathlib
import pandas as pd
import os
import time

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
configfile: "config/config.yaml"

# setup
run_path = config["run"]
run_name = run_path.split("/")[-1]

archive_path = config["archive"]

bcl2fastq_out_path = "{}/Run/{}/{}".format(archive_path, run_name, config["bcl2fastq"])
log_path = "{}/Log/{}".format(archive_path, config["log"])
os.makedirs(log_path, exist_ok = True)

samplesheet_path = config["samplesheet"]
sample_info = extract_sample_info(samplesheet_path)
samples = list(sample_info.loc[:,"Sample_ID"])
sample_info.index = samples

################################################################################

rule all:
  input:
    bcl2fastq_out_path + "/Stats/Stats.json"

rule bcl2fastq:
  input:
    samplesheet_path,
    run_path + "/RTAComplete.txt"
  output:
    bcl2fastq_out_path + "/Stats/Stats.json"
  params:
    name = "bcl2fastq",
    thread = "6",
    mem_per_thread = "6G",
    output = "{}/bcl2fastq.output".format(log_path),
    error = "{}/bcl2fastq.error".format(log_path)
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
