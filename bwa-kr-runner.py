# This script runs BWA-MEM2, Samtools view&fastq, Kraken2 and taxpasta on files from a samplesheet
# Make sure you are in an active conda environment with bwa-mem2, samtools, kraken2 and taxpasta installed

import csv
import subprocess
import os

# Path to bwa-mem2 index of the reference genome
    # Modify if needed
bwamem2_db = "/home/birgit/data/db_index/bwa-mem2/human"

# Path to kraken2 database
    # Modify if needed
kraken2_db = "/home/birgit/data/db_index/kraken2/standard"


# Add path to taxpasta taxonomy (kraken2) database!!

# Path to CSV file containing sample information
    # Modify if neeeded 
samplesheet = "/home/birgit/data/testfiles/liu_metagenomics/liu_meta_samplesheet_resume.csv"

# Path to dir where the output files will be stored
    # Modify if needed
output_dir = "/mnt/FS1/data_1/PRJNA866654_Lung_Cancer_Gut_Microbiome/Results/"

# Create output dir if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Create subdirectories for bwa-mem2 and kraken2 output if they do not exist
if not os.path.exists(f"{output_dir}/bwamem2"):
    os.makedirs(f"{output_dir}/bwamem2")
if not os.path.exists(f"{output_dir}/kraken2"):
    os.makedirs(f"{output_dir}/kraken2")

# Function to run BWA-MEM2 for each sample
def run_bwa_mem2(sample, fastq_1, fastq_2):
    # Construct the command to run BWA-MEM2
    command = [
        "bwa-mem2",
        "mem",
        "-t",
        "32",
        bwamem2_db,
        fastq_1,
        fastq_2,
    ]
    # Create file for output
    output_file = f"{output_dir}/bwamem2/{sample}_out.sam"
    # Open output file
    with open(output_file, "w") as f_out:
        # Run the BWA-MEM2 command and write output to file
        subprocess.run(command, stdout=f_out)
    
    # Print status message
    print(f"\nFinished BWA-MEM2 for sample {sample}\n")

# Function to run samtools for each sample
def run_samtools(sample):
    # Construct the command to run samtools view
    command_view = [
        "samtools",
        "view",
        "-h",
        "-F2",
        f"{output_dir}/bwamem2/{sample}_out.sam"
    ]
    # Create file for output
    output_file = f"{output_dir}/bwamem2/{sample}_out_unconc.sam"
    # Open output file
    with open(output_file, "w") as f_out:
        # Run the samtools view command and write output to file
        subprocess.run(command_view, stdout=f_out)

    # Construct the command to run samtools fastq
    command_fastq = [
        "samtools",
        "fastq",
        f"{output_dir}/bwamem2/{sample}_out_unconc.sam",
        "-1",
        f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq",
        "-2",
        f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq"
    ]
    # Run the samtools fastq command
    subprocess.run(command_fastq)

    # Print status message
    print(f"\nFinished samtools view and fastq for sample {sample}\n")

def run_kraken2(sample):
    # Command to run kraken2s
    command = [
        "kraken2",
        "--threads",
        "32",
        "--db",
        kraken2_db,
        "--paired",
        f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq",
        f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq",
        "--output",
        f"{output_dir}/kraken2/{sample}_kraken2.out",
        "--report",
        f"{output_dir}/kraken2/{sample}_kraken2_report.txt"
    ]
    # Run the kraken2 command
    subprocess.run(command)
    
    # Print status message
    print(f"\nFinished kraken2 for sample {sample}\n")

def run_taxpasta():
    # Command to run taxpasta
    command = [
        "taxpasta",
        "merge",
        "--profiler",
        "kraken2",
        "--output",
        f"{output_dir}/merged_samples.biom",
        "--output-format",
        "BIOM",
        "--taxonomy",
        "/home/birgit/data/db_index/kraken2/standard/taxonomy/"
        "--samplesheet",
        f"{output_dir}/samplesheet_kraken_reports.csv"
    ]
    # Run the taxpasta command
    subprocess.run(command)
    
    # Print status message
    print("\nFinished taxpasta\n")

# Intitalize list to store the sample names and kraken2 report locations
samplesheet_kraken_reports = []

# Read sample information from the CSV file and run tools for each sample
with open(samplesheet, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row["sample"]
        fastq_1 = row["fastq_1"]
        fastq_2 = row["fastq_2"]
        # Run BWA-MEM2 for the current sample
        run_bwa_mem2(sample, fastq_1, fastq_2)
        # Run samtools for the current sample
        run_samtools(sample)
        # Run kraken2 for the current sample
        run_kraken2(sample)
        # Append the sample name and kraken2 report location to the list
        samplesheet_kraken_reports.append({
            "sample": sample,
            "profile": f"{output_dir}/kraken2/{sample}_kraken2_report.txt"
        })

# Write the sample sheet to a CSV file
with open(f"{output_dir}/samplesheet_kraken_reports.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["sample", "profile"])
    writer.writeheader()
    writer.writerows(samplesheet_kraken_reports)

# Run taxpasta to merge the kraken2 reports
run_taxpasta()

# Last status message
print(f"\nFinished all {len(samplesheet_kraken_reports)} samples :)\n")
