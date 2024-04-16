"""
This script counts reads in sam outputs and kraken2 reports.

Make sure you are in an active conda environment with samtools installed.

Usage:
python counter.py --samplesheet <path_to_samplesheet> --out_dir <output_directory>

Arguments:
--samplesheet: Path to the CSV file containing sample information.
--out_dir: Path to the output directory (not ending with /).

The script performs the following steps for each sample in the samplesheet:
1.


Example usage:
python counter.py --samplesheet /path/to/samplesheet.csv --out_dir /path/to/output_directory
"""

# Import modules
import csv
import subprocess
import os
import argparse
import re


def flagstat(sample, output_dir, paired):
    """
    Runs samtools flagstat command on a given sample, extracts number of input reads, concordantly and unconcordantly aligned reads.
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
    if paired:
        # Select the 12th line (amount of properly paired reads)
        twelveth_line = lines[11]
        # Extract the number from the 12th line
        properly_paired = int(twelveth_line.split(' ')[0])
        # Calculate amount of disconcordantly paired reads
        unconcordantly_aligned = total_input - properly_paired
    else:
        # Select the 8th line (amount of properly paired reads)
        eight_line = lines[7]
        # Extract the number from the 12th line
        properly_paired = int(eight_line.split(' ')[0])
        # Calculate amount of disconcordantly paired reads
        unconcordantly_aligned = total_input - properly_paired
    # Print the results
    print(f"Number of BWA-MEM2 input reads for sample {sample}: {total_input}", flush=True)
    print(f"Number of concordantly aligned reads for sample {sample}: {properly_paired}", flush=True)
    print(f"Number of unaligned reads for sample {sample}: {unconcordantly_aligned}", flush=True)
    # Return amounts
    return total_input, properly_paired, unconcordantly_aligned

# function to extract amount of homo sapiens and bacteria kraken2 classifications from kreport
def kraken2_stats(sample, output_dir):
    """
    Calculate the number of Homo sapiens and Bacterial classified reads for a given sample based on a Kraken2 report.
    Prints status message to terminal with the results, warns when no classification count is zero.

    Args:
        sample (str): The name of the sample.
        output_dir (str): The output directory where the Kraken2 report is located.

    Returns:
        tuple: A tuple containing the number of Homo sapiens reads and the number of Bacterial reads.

    Raises:
        FileNotFoundError: If the Kraken2 report file is not found.
    """
    # Construct filename for report of sample
    filename = f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"

    # Intitalize counters (outside of for loop, so stays 0 when nothing is encountered and no errors)
    homo_number = 0
    bac_number = 0
    input_number = 0
    unclassfied_number = 0
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
            # If Bacteria is encountered (based on tax id, some species contain 'Bacteria' in name), take value and print
            if re.fullmatch('^2$', columns[4]):
                bac_number = int(columns[1])
                print(f"Number of Bacterial reads for sample {sample}: {bac_number}", flush=True)
            if re.fullmatch('^0$', columns[4]):
                unclassfied_number = int(columns[1])
                print(f"Number of unclassified reads for sample {sample}: {unclassfied_number}", flush=True)
    # Warn user if no Homo sapiens classifications
    if homo_number == 0:
        print(f"Warning: No Homo sapiens classifications found in kraken2 report of sample {sample}.\n", flush=True)
    # Warn user if no bacteria classifications found
    if bac_number == 0:
        print(f"Warning: No Bacteria classifications found in kraken2 report of sample {sample}.\n", flush=True)
    # Return values
    return homo_number, bac_number, input_number, unclassfied_number

# Function that creates csv based on outputs of flagstat and extract homo sapiens reads functions
def read_count_extractor(sample, output_dir, paired):
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
    output_file = f"{output_dir}/read_counts.csv"
    # Run flagstat for the current sample
    total_input, properly_paired, unconcordantly_aligned = flagstat(sample, output_dir, paired)
    # extract homo sapiens reads
    homo_number, bac_number, input_number, unclassifed_number = kraken2_stats(sample, output_dir)
    # Write values from each sample to csv
    with open(output_file, 'a', newline='') as outfile:
        out_writer = csv.writer(outfile)
        out_writer.writerow([sample, total_input, properly_paired, unconcordantly_aligned, input_number, homo_number, bac_number, unclassifed_number])

# Main function that calls other functions based on samplesheet and passed arguments,
def main():
    """
    """
    # Create parser
    parser = argparse.ArgumentParser(description='Run counter.')
    # Add the arguments
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to CSV file containing sample information')
    parser.add_argument('--out_dir', type=str, required=True, help='Path to output directory (not ending with /)')
    # Parse!
    args = parser.parse_args()

    # Assign arguments to variables
    samplesheet = args.samplesheet
    output_dir = args.out_dir

    # Create file for summary of extracted readcounts
    output_file = f"{output_dir}/read_counts.csv"

    # Create output file if it does not exist
    with open(output_file, 'w', newline='') as outfile:
        # Create headers
        headers = ['sample', 'total_input_reads', 'human_bwa', 'not_human_bwa', 'input_kraken', 'human_kraken', 'bacterial_kraken', 'unclassified_kraken']
        # Write headers to the output file
        out_writer = csv.writer(outfile)
        out_writer.writerow(headers)
    
    # Initialize a counter for the number of samples processed
    num_samples_processed = 0

    # Read sample information from the CSV file and call functions for each sample
    with open(samplesheet, "r") as ssheet:
        ssheet_reader = csv.DictReader(ssheet)
        paired = "fastq_2" in ssheet_reader.fieldnames
        for row in ssheet_reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"] if paired else None
            # Extract read counts and put in summary csv
            read_count_extractor(sample, output_dir, paired)
            # Up the counter
            num_samples_processed += 1

    print(f"\nCreated {output_file} with summarized read counts of all {num_samples_processed} samples :)\n", flush=True)

# If this script is run directly call the main function
if __name__ == '__main__':
    main()
