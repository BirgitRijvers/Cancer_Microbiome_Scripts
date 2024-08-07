"""
This script runs Kraken2 on files from a samplesheet.
It detects if the input data is single or paired end based on the samplesheet columns. 

Make sure you are in an active conda environment with Kraken2 installed.

Usage:
python kr_runner.py --samplesheet <path_to_samplesheet> --out_dir <output_directory>

Arguments:
--samplesheet: Path to the CSV file containing sample information.
--out_dir: Path to the output directory (not ending with /).

The script performs the following steps for each sample in the samplesheet:
1. Runs Kraken2 to classify the unaligned FASTQ files and generate a report. (path to database hardcoded in script)

The script also creates subdirectories for Kraken2 output if they do not exist.
The output directory should not end with a forward slash (/).

Note: Modify the paths to the Kraken2 database if needed.

Example usage:
python pipeline.py --samplesheet /path/to/samplesheet.csv --out_dir /path/to/output_directory
"""

# Import modules
import csv
import subprocess
import os
import argparse

# Path to kraken2 database
    # Modify if neeeded
kraken2_db = "/mnt/FS2/data_2/Users/Birgit/db_index/kraken2/standard/"

def run_kraken2(sample, fastq_1, fastq_2, output_dir):
    """
    Runs Kraken2 with 64 threads for a given sample.
    Creates Kraken2 output file and report in reports directory
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
        "64",
        "--db",
        kraken2_db,
        "--paired",
        fastq_1,
        fastq_2,
        "--output",
        f"{output_dir}/kraken2/{sample}_kraken2.out",
        "--report",
        f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"
    ]

    subprocess.run(command_kraken2_paired)
    # # Check if paired or single end data
    # if os.path.exists(f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"):
    #     # Run the kraken2 command for paired end data
    #     subprocess.run(command_kraken2_paired)
    # else:
    #     # Run the kraken2 command for single end data
    #     subprocess.run(command_kraken2_single)

    # Print status message
    print(f"\nFinished kraken2 for sample {sample}\n", flush=True)

# Main function that calls other functions based on samplesheet and passed arguments,
def main():
    """
    Main function that calls other functions based on samplesheet and passed arguments.
    """
    # Create parser
    parser = argparse.ArgumentParser(
        description=(
            'Runs Kraken2 on all provided samples.'
        )
    )
    # Add the arguments
    parser.add_argument('--samplesheet', type=str, required=True,
                        help='Path to CSV file containing sample information')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='Path to Kraken2 output directory (not ending with /)')
    # Parse!
    args = parser.parse_args()

    # Assign arguments to variables
    samplesheet = args.samplesheet
    output_dir = args.out_dir

    # Create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create subdirectories for Kraken2 output
    if not os.path.exists(f"{output_dir}/kraken2"):
        os.makedirs(f"{output_dir}/kraken2")
    if not os.path.exists(f"{output_dir}/kraken2/reports"):
        os.makedirs(f"{output_dir}/kraken2/reports")

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
            run_kraken2(sample, fastq_1, fastq_2, output_dir)
            num_samples_processed += 1

    # Status messages
    print(f"\nProcessed all {num_samples_processed} samples :)\n", flush=True)

# If this script is run directly call the main function
if __name__ == '__main__':
    main()
