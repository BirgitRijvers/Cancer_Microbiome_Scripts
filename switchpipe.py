"""
This script runs fastp, BWA-MEM2, Samtools view&fastq, Kraken2 and kraken-biom on files from a samplesheet.
It also creates a CSV file with read counts from the SAM files and Kraken reports.
It detects if the input data is single or paired end based on the samplesheet columns. 

Make sure you are in an active conda environment with fastp, bwa-mem2, samtools, kraken2, and kraken-biom installed.

Usage:
python pipeline.py --samplesheet <path_to_samplesheet> --out_dir <output_directory>

Arguments:
--samplesheet: Path to the CSV file containing sample information.
--out_dir: Path to the output directory (not ending with /).

The script performs the following steps for each sample in the samplesheet:
1. Runs fastp for QC, trimming and filtering of the FASTQ files.
2. Runs BWA-MEM2 to align the QCed FASTQ files to a reference genome. (path to database hardcoded in script)
3. Runs Samtools view to filter and extract unaligned reads from the BWA-MEM2 output.
4. Runs Samtools fastq to convert the filtered SAM file to FASTQ files.
5. Runs Kraken2 to classify the unaligned FASTQ files and generate a report. (path to database hardcoded in script)
6. Extracts the number of Homo sapiens and bacterial classifications from the Kraken2 report.
7. Runs Samtools flagstat to calculate read statistics from the BWA-MEM2 output.
8. Writes the extracted read counts from BWA-MEM2 and Kraken2 to a summary CSV file.
9. Runs Kraken-biom to convert the Kraken2 reports to a BIOM file.

The script also creates subdirectories for BWA-MEM2 and Kraken2 output if they do not exist.
The output directory should not end with a forward slash (/).

Note: Modify the paths to the BWA-MEM2 index and Kraken2 database if needed.

Example usage:
python pipeline.py --samplesheet /path/to/samplesheet.csv --out_dir /path/to/output_directory
"""

# Import modules
import csv
import subprocess
import os
import argparse
import re
import json

# Path to bwa-mem2 index of the reference genome
    # Modify if needed
bwamem2_db = "/mnt/FS2/data_2/Users/Birgit/db_index/bwa-mem2/human"

# Path to kraken2 database
    # Modify if neeeded
kraken2_db = "/mnt/FS2/data_2/Users/Birgit/db_index/kraken2/standard/"


def run_fastp(sample, fastq_1, fastq_2, output_dir):
    """
    Run fastp command for quality control, trimming and filtering of sequencing data.

    Args:
        sample (str): The sample name.
        fastq_1 (str): Path to the first input FASTQ file.
        fastq_2 (str): Path to the second input FASTQ file, None if single end data.
        output_dir (str): Path to the output directory.

    Returns:
        tuple or str: If `fastq_2` is provided, returns a tuple containing the paths to the output FASTQ files for paired data.
                      If `fastq_2` is None, returns the path to the output FASTQ file for single end data.
    """
    # Paths to output fastq files for paired data
    qc_fastq_1 = f"{output_dir}/fastp/{sample}_1.fastq"
    qc_fastq_2 = f"{output_dir}/fastp/{sample}_2.fastq"

    # Construct paired end command
    command_paired = [
        "fastp",
        "-i",
        fastq_1,
        "-I",
        fastq_2,
        "-o",
        f"{output_dir}/fastp/{sample}_1.fastq",
        "-O",
        f"{output_dir}/fastp/{sample}_2.fastq",
        "-h",
        f"{output_dir}/fastp/reports/html/{sample}.html",
        "-j",
        f"{output_dir}/fastp/reports/json/{sample}.json",
        "-R",
        f"fastp report {sample}",
        "--detect_adapter_for_pe"
    ]

    # Path to the output fastq file single end data
    qc_fastq = f"{output_dir}/fastp/{sample}.fastq"
    
    # Construct single end command
    command_single = [
        "fastp",
        "-i",
        fastq_1,
        "-o",
        f"{output_dir}/fastp/{sample}.fastq",
        "-h",
        f"{output_dir}/fastp/reports/html/{sample}.html",
        "-j",
        f"{output_dir}/fastp/reports/json/{sample}.json",
        "-R",
        f"fastp report {sample}"
    ]
    
    # Check if paired or single end data
    if fastq_2:
        # Run fastp command for paired end data
        subprocess.run(command_paired)
        return qc_fastq_1, qc_fastq_2
    else:
        # Run fastp command for single end data
        subprocess.run(command_single)
        # set qc_fastq_2 to None for single end data
        qc_fastq_2 = None
        return qc_fastq, qc_fastq_2

def run_bwa_mem2(sample, fastq_1, fastq_2, output_dir):
    """
    Run BWA-MEM2 for a given sample using the provided FASTQ filess.
    Creates sam file with BWA-MEM2 output.
    Print status message to terminal when finished.

    Args:
        sample (str): The name of the sample.
        fastq_1 (str): The path to the first paired FASTQ file.
        fastq_2 (str): The path to the second paired FASTQ file, None for single end data.
        output_dir (str): The directory where the output sam file will be saved.

    Returns:
        None
    """

    # Construct the command to run BWA-MEM2 on paired end data
    command_paired = [
        "bwa-mem2",
        "mem",
        "-t",
        "32",
        bwamem2_db,
        fastq_1,
        fastq_2,
    ]
    # Construct the command to run BWA-MEM2 on single end data
    command_single = [
        "bwa-mem2",
        "mem",
        "-t",
        "32",
        bwamem2_db,
        fastq_1
    ]
    # Create file for output
    output_file = f"{output_dir}/bwamem2/{sample}_out.sam"
    # Open output file
    with open(output_file, "w") as f_out:
        if fastq_2:
            # Run the BWA-MEM2 command for paired end data and write output to file
            subprocess.run(command_paired, stdout=f_out)
        else:
            # Run the BWA-MEM2 command for single end data and write output to file
            subprocess.run(command_single, stdout=f_out)
    # Print status message
    print(f"\nFinished BWA-MEM2 for sample {sample}\n", flush=True)

def run_samtools(sample, output_dir, fastq_2):
    """
    Runs samtools view and samtools fastq commands on BWA-MEM2 output on given sample.
    Creates sam file with unaligned reads and fastq files with unaligned reads.
    Print status message to terminal when finished with each command.

    Args:
        sample (str): The name of the sample, from samplesheet.
        output_dir (str): The output directory where the sam and fastq files will be saved.

    Returns:
        None
    """
    # Construct the commands to run samtools view
    # -f12 to only extract unaligned reads of paired end data
    command_view_paired = [
        "samtools",
        "view",
        "-h",
        "-f12",
        f"{output_dir}/bwamem2/{sample}_out.sam"
    ]
    
    # -f4 to only extract unaligned reads of single end data
    command_view_single = [
    "samtools",
    "view",
    "-h",
    "-f4",
    f"{output_dir}/bwamem2/{sample}_out.sam"
    ]
    
    # Create file for output
    output_file = f"{output_dir}/bwamem2/{sample}_out_unconc.sam"
    # Open output file
    with open(output_file, "w") as f_out:
        # Check if paired or single end data
        if fastq_2:
            # Run the samtools view command for paired data and write output to file
            subprocess.run(command_view_paired, stdout=f_out)
        else:
            # Run the samtools view command for single data and write output to file
            subprocess.run(command_view_single, stdout=f_out)
    
    # Print status message
    print(f"\nFinished samtools view for sample {sample}\n", flush=True)

    # Construct the command to run samtools fastq paired end data
    command_fastq_paired = [
        "samtools",
        "fastq",
        f"{output_dir}/bwamem2/{sample}_out_unconc.sam",
        "-1",
        f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq",
        "-2",
        f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq"
    ]
    
    # Construct the command to run samtools fastq single end data
    command_fastq_single = [
        "samtools",
        "fastq",
        f"{output_dir}/bwamem2/{sample}_out_unconc.sam"
    ]
    
    # Run the samtools fastq command based on paired or single end data
    if fastq_2:
        # Run the samtools fastq command for paired end data and write output to file
        subprocess.run(command_fastq_paired)
    else:
        # Create file for output single end data
        output_fastq = f"{output_dir}/bwamem2/{sample}_out_unconc.fastq"
        # Open output file
        with open(output_fastq, "w") as fastq_out:
            # Run the samtools fastq command for single end data and write output to file
            subprocess.run(command_fastq_single, stdout=fastq_out)

    # Print status message
    print(f"\nFinished samtools fastq for sample {sample}\n", flush=True)

def check_fastq(sample, output_dir):
    """
    Checks if the fastq files contain data, exits with error if it is empty.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory.

    Returns:
        None
    
    Raises:
        ValueError: If the fastq file is empty
    """
    # Check if paired or single data
    if os.path.exists(f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"):
        # Construct filenames for fastq files
        fastq_file_1 = f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"
        fastq_file_2 = f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq"
        # Check if fastqfiles are empty
        if os.path.getsize(fastq_file_1) == 0:
            raise ValueError(f"Error: Fastq file {fastq_file_1} is empty for sample {sample}.\n")
        if os.path.getsize(fastq_file_2) == 0:
            raise ValueError(f"Error: Fastq file {fastq_file_2} is empty for sample {sample}.\n")
    else:
        fastq_file_single = f"{output_dir}/bwamem2/{sample}_out_unconc.fastq"
        # Check if fastqfile is empty
        if os.path.getsize(fastq_file_single) == 0:
            raise ValueError(f"Error: Fastq file {fastq_file_single} is empty for sample {sample}.\n")

def run_kraken2(sample, output_dir):
    """
    Runs Kraken2 with 32 threads for a given sample.
    Creates kraken2 output file and report in reports directory
    Print status message to terminal when finished.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory.

    Returns:
        None
    """
    # Command to run kraken2 on paired end data
    command_kraken2_paired = [
        "kraken2",
        "--threads",
        "32",
        "--db",
        kraken2_db,
        "--paired",
        f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq",
        f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq",
        "--confidence",
        "0.05",
        "--output",
        f"{output_dir}/kraken2/{sample}_kraken2.out",
        "--report",
        f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"
    ]

    # Command to run kraken2 on single end data
    command_kraken2_single = [
        "kraken2",
        "--threads",
        "32",
        "--db",
        kraken2_db,
        f"{output_dir}/bwamem2/{sample}_out_unconc.fastq",
        "--output",
        f"{output_dir}/kraken2/{sample}_kraken2.out",
        "--report",
        f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"
    ]

    # Check if paired or single end data
    if os.path.exists(f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"):
        # Run the kraken2 command for paired end data
        subprocess.run(command_kraken2_paired)
    else:
        # Run the kraken2 command for single end data
        subprocess.run(command_kraken2_single)

    # Print status message
    print(f"\nFinished kraken2 for sample {sample}\n", flush=True)

def fastp_extraction(sample, output_dir):
    """
    Extracts the total number of reads from a JSON report file generated by fastp.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory where the JSON report file is located.

    Returns:
        int: The total number of reads for the given sample.
    """
    # Construct filename for JSON report of sample
    filename = f"{output_dir}/fastp/reports/json/{sample}.json"

    # Open JSON file and load it into a dictionary
    with open(filename, 'r') as file:
        report = json.load(file)
    # Extract total number of reads
    total_reads = report['summary']['before_filtering']['total_reads']

    # Print status message
    print(f"\nTotal number of input reads for sample {sample}: {total_reads}", flush=True)

    # Return total number of reads
    return total_reads

def flagstat(sample, output_dir):
    """
    Runs samtools flagstat command on a given sample, extracts number of input reads, 
    concordantly and unconcordantly aligned reads.
    Prints status message to terminal with the results.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory.

    Returns:
        tuple: A tuple containing the total number of input reads, 
               the number of concordantly aligned reads,
               and the number of unconcordantly aligned reads.
    """
    # Construct command
    command = [
        "samtools",
        "flagstat",
        "-@",
        "32",
        f"{output_dir}/bwamem2/{sample}_out.sam"
    ]

    # Run samtools flagstat command
    result = subprocess.run(command, capture_output=True, text=True)
    # Split the output into lines
    lines = result.stdout.split('\n')

    # Select the second line (amount of input reads)
    second_line = lines[1]
    # Extract the amount of input reads from 2nd line
    total_input = int(second_line.split(' ')[0])

    # Select the 12th line (amount of properly paired end reads)
    twelveth_line = lines[11]
    # Extract the number from the 12th line
    properly_mapped_paired = int(twelveth_line.split(' ')[0])

    # Select the 8th line (primary mapped single end reads)
    eight_line = lines[7]
    # Extract the number from the 8th line
    properly_mapped_single = int(eight_line.split(' ')[0])

    # Print the amount of input reads
    print(f"Number of BWA-MEM2 input reads for sample {sample}: {total_input}", flush=True)
    
    # Calculate amount of disconcordantly paired reads
    # Check if paired end data
    if os.path.exists(f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"):
         # Calculate amount of unaligned paired reads
        unaligned_paired = total_input - properly_mapped_paired
        # Print results
        print(f"Number of concordantly aligned reads for sample {sample}: {properly_mapped_paired}", flush=True)
        print(f"Number of unaligned reads for sample {sample}: {unaligned_paired}", flush=True)
        # Return amounts
        return total_input, properly_mapped_paired, unaligned_paired
    # Single end data
    else:
        # Calculate amount of unaligned single reads
        unaligned_single = total_input - properly_mapped_single
        # Print results
        print(f"Number of concordantly aligned reads for sample {sample}: {properly_mapped_single}", flush=True)
        print(f"Number of unaligned reads for sample {sample}: {unaligned_single}", flush=True)
        # Return amounts
        return total_input, properly_mapped_single, unaligned_single

def kraken2_stats(sample, output_dir):
    """
    Calculate classifications of reads for a given sample based on a Kraken2 report.
    Prints status message to terminal with the results, 
    warns when human/bacterial classification count is zero.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory where the Kraken2 report is located.

    Returns:
        tuple: A tuple containing the number of input Homo sapiens, Bacterial reads 
        and unclassified reads.

    Raises:
        FileNotFoundError: If the Kraken2 report file is not found.
    """
    # Construct filename for report of sample
    filename = f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"

    # Intitalize counters (outside of loop, so stays 0 when nothing is encountered and no errors)
    homo_number = 0
    bac_number = 0
    unclassified_number = 0
    input_number = 0
    
    # Open report, loop through lines
    with open(filename, 'r') as file:
        for line in file:
            # Split columns based on tab
            columns = line.split('\t')
            # Calculate total number
            input_number += int(columns[2])
            # If homo sapiens is encountered, take value and print
            if 'Homo sapiens' in columns[5]:
                homo_number = int(columns[1])
                print(f"Number of Homo sapiens reads for sample {sample}: {homo_number}", flush=True)
            # If Bacteria is encountered (based on tax id, some species contain 'Bacteria' in name), take value, print
            if re.fullmatch('^2$', columns[4]):
                bac_number = int(columns[1])
                print(f"Number of Bacterial reads for sample {sample}: {bac_number}", flush=True)
            # If Unclassified is found (based on tax id, some species contain 'Unclassified' in name), take value, print
            if re.fullmatch('^0$', columns[4]):
                unclassified_number = int(columns[1])
                print(f"Number of unclassified reads for sample {sample}: {unclassified_number}", flush=True)            
    
    # Warn user if no Homo sapiens classifications are found
    if homo_number == 0:
        print(f"Warning: No Homo sapiens classifications found in kraken2 report of sample {sample}.\n", flush=True)
    
    # Warn user if no Bacteria classifications found
    if bac_number == 0:
        print(f"Warning: No Bacteria classifications found in kraken2 report of sample {sample}.\n", flush=True)
    
    # Return values
    return homo_number, bac_number, unclassified_number, input_number

def read_count_extractor(sample, output_dir):
    """
    Extracts read counts for a given sample and writes them to a CSV file.
    Calls kraken2_stats and flagstat functions.
    Writes information to summary CSV file ("read_counts.csv").

    Args:
        sample (str): The path to the input sample file.
        output_dir (str): The directory where the output CSV file will be written.

    Returns:
        None
    """
    # Create path to output file
    output_file = f"{output_dir}/read_counts.csv"
    
    # Run fastp_summary for the current sample
    raw_input = fastp_extraction(sample, output_dir)
    
    # Run flagstat for the current sample
    bwa_input, properly_paired, unconcordantly_aligned = flagstat(sample, output_dir)
    
    # Extract read counts from kraken2 report with kraken2_stats
    homo_number, bac_number, unclassified_number, total_number = kraken2_stats(sample, output_dir)
    
    # Write values from each sample to csv
    with open(output_file, 'a', newline='') as outfile:
        out_writer = csv.writer(outfile)
        out_writer.writerow([sample, raw_input, 
                             bwa_input, properly_paired, unconcordantly_aligned, 
                             total_number, homo_number, bac_number, unclassified_number])

def run_kraken_biom(output_dir):
    """
    Runs kraken-biom command to generate a classified.biom file with taxonomic information of all samples.
    Prints status message to terminal when finished.

    Parameters:
    output_dir (str): The output directory where the classified.biom file will be saved.

    Returns:
    None
    """
    # Command to run kraken-biom
    command_kr_biom = [
        "kraken-biom",
        "-k",
        f"{output_dir}/kraken2/reports/",
        "-o",
        f"{output_dir}/classified.biom"
    ]
    
    # Run kraken-biom command
    subprocess.run(command_kr_biom)
    # Print status message
    print(f"\nCreated {os.path.join(output_dir,'classified.biom')} with taxonomic information of all samples :)\n", flush=True)

# Main function that calls other functions based on samplesheet and passed arguments,
def main():
    """
    Run pipeline with BWA-MEM2 host depletion, samtools, kraken2, 
    biom conversion and read count summarizer.
    Does this by calling other functions.

    Args:
        --samplesheet (str): Path to CSV file containing sample information.
        --out_dir (str): Path to output directory (not ending with /).

    Returns:
        None
    """
    # Create parser
    parser = argparse.ArgumentParser(
        description=(
            'Run pipeline with BWA-MEM2 host depletion, samtools, kraken2, '
            'biom conversion and read count summarizer.'
        )
    )
    # Add the arguments
    parser.add_argument('--samplesheet', type=str, required=True,
                        help='Path to CSV file containing sample information')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='Path to output directory (not ending with /)')
    # Parse!
    args = parser.parse_args()

    # Assign arguments to variables
    samplesheet = args.samplesheet
    output_dir = args.out_dir

    # Create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create subdirectories for qc, bwa-mem2 and kraken2 output if they do not exist
    if not os.path.exists(f"{output_dir}/fastp"):
        os.makedirs(f"{output_dir}/fastp")
    if not os.path.exists(f"{output_dir}/fastp/reports"):
        os.makedirs(f"{output_dir}/fastp/reports")
    if not os.path.exists(f"{output_dir}/fastp/reports/html"):
        os.makedirs(f"{output_dir}/fastp/reports/html")
    if not os.path.exists(f"{output_dir}/fastp/reports/json"):
        os.makedirs(f"{output_dir}/fastp/reports/json")
    if not os.path.exists(f"{output_dir}/bwamem2"):
        os.makedirs(f"{output_dir}/bwamem2")
    if not os.path.exists(f"{output_dir}/kraken2"):
        os.makedirs(f"{output_dir}/kraken2")
    if not os.path.exists(f"{output_dir}/kraken2/reports"):
        os.makedirs(f"{output_dir}/kraken2/reports")

    # Create file for summary of extracted readcounts
    readcount_file = f"{output_dir}/read_counts.csv"

    # Create output file if it does not exist
    with open(readcount_file, 'w', newline='') as outfile:
        # Create headers
        headers = ['sample', 'total_input_reads',
                   'bwa_input', 'human_bwa', 'not_human_bwa',
                   'input_kraken', 'human_kraken', 'bacterial_kraken', 'unclassified_kraken']
        # Write headers to the output file
        out_writer = csv.writer(outfile)
        out_writer.writerow(headers)

    # Initialize a counter for the number of samples processed
    num_samples_processed = 0

    # Read sample information from the CSV file and call functions for each sample
    with open(samplesheet, "r") as ssheet:
        ssheet_reader = csv.DictReader(ssheet)
        # Determine if data is paired or single end
        paired = "fastq_2" in ssheet_reader.fieldnames
        for row in ssheet_reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            # Set fastq_2 to None if single end data	
            fastq_2 = row["fastq_2"] if paired else None
            # Run fastp for the current sample
            fastq_1, fastq_2 = run_fastp(sample, fastq_1, fastq_2, output_dir)
            # Run BWA-MEM2 for the current sample
            run_bwa_mem2(sample, fastq_1, fastq_2, output_dir)
            # Run samtools for the current sample
            run_samtools(sample, output_dir, fastq_2)
            # Check fastq file
            check_fastq(sample, output_dir)
            # Run kraken2 for the current sample
            run_kraken2(sample, output_dir)
            # Extract read counts and put in summary csv
            read_count_extractor(sample, output_dir)
            # Up the counter
            num_samples_processed += 1

    # Run kraken-biom
    run_kraken_biom(output_dir)

    # Status messages
    print(f"\nProcessed all {num_samples_processed} samples :)\n", flush=True)
    print(f"\nCreated {readcount_file} with summarized read counts of all {num_samples_processed} samples :)\n", flush=True)

# If this script is run directly call the main function
if __name__ == '__main__':
    main()
