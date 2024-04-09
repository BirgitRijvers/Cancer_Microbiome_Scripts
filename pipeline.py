"""
This script runs BWA-MEM2, Samtools view&fastq, Kraken2 and kraken-biom on files from a samplesheet.
It also creates a CSV file with read counts from the SAM files and Kraken reports. 

Make sure you are in an active conda environment with bwa-mem2, samtools, kraken2, and kraken-biom installed.

Usage:
python pipeline.py --samplesheet <path_to_samplesheet> --out_dir <output_directory>

Arguments:
--samplesheet: Path to the CSV file containing sample information.
--out_dir: Path to the output directory (not ending with /).

The script performs the following steps for each sample in the samplesheet:
1. Runs BWA-MEM2 to align the paired-end FASTQ files to a reference genome. (path to database hardcoded in script)
2. Runs Samtools view to filter and extract unaligned reads from the BWA-MEM2 output.
3. Runs Samtools fastq to convert the filtered SAM file to FASTQ files.
4. Runs Kraken2 to classify the unaligned FASTQ files and generate a report. (path to database hardcoded in script)
5. Extracts the number of Homo sapiens and bacterial classifications from the Kraken2 report.
6. Runs Samtools flagstat to calculate read statistics from the BWA-MEM2 output.
7. Writes the extracted read counts from BWA-MEM2 and Kraken2 to a summary CSV file.
8. Runs Kraken-biom to convert the Kraken2 reports to a BIOM file.

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

# Path to bwa-mem2 index of the reference genome
    # Modify if needed
bwamem2_db = "/home/birgit/data/db_index/bwa-mem2/human"

# Path to kraken2 database
    # Modify if neeeded
kraken2_db = "/home/birgit/data/db_index/kraken2/standard"

def run_bwa_mem2(sample, fastq_1, fastq_2, output_dir):
    """
    Run BWA-MEM2 for a given sample using the provided FASTQ files.
    Creates sam file with BWA-MEM2 output.
    Print status message to terminal when finished.

    Args:
        sample (str): The name of the sample.
        fastq_1 (str): The path to the first paired FASTQ file.
        fastq_2 (str): The path to the second paired FASTQ file.
        output_dir (str): The directory where the output sam file will be saved.

    Returns:
        None
    """
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
    print(f"\nFinished BWA-MEM2 for sample {sample}\n", flush=True)

def run_samtools(sample, output_dir):
    """
    Runs samtools view -f12 and samtools fastq commands on BWA-MEM2 output on given sample.
    Creates sam file with unaligned reads and fastq files with unaligned reads.
    Print status message to terminal when finished with each command.

    Args:
        sample (str): The name of the sample, from samplesheet.
        output_dir (str): The output directory where the sam and fastq files will be saved.

    Returns:
        None
    """
    # Construct the command to run samtools view
    # Use -f12 to only extract unaligned reads
    command_view = [
        "samtools",
        "view",
        "-h",
        "-f12",
        f"{output_dir}/bwamem2/{sample}_out.sam"
    ]
    # Create file for output
    output_file = f"{output_dir}/bwamem2/{sample}_out_unconc.sam"
    # Open output file
    with open(output_file, "w") as f_out:
        # Run the samtools view command and write output to file
        subprocess.run(command_view, stdout=f_out)
    # Print status message
    print(f"\nFinished samtools view for sample {sample}\n", flush=True)

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
    subprocess.run(command_fastq)
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
    # Path to fastq file
    fastq_file_1 = f"{output_dir}/bwamem2/{sample}_out_unconc_1.fastq"
    fastq_file_2 = f"{output_dir}/bwamem2/{sample}_out_unconc_2.fastq"
    # Check if fastqfile is empty
    if os.path.getsize(fastq_file_1) == 0:
        raise ValueError(f"Error: Fastq file {fastq_file_1} is empty for sample {sample}.\n")
    if os.path.getsize(fastq_file_2) == 0:
        raise ValueError(f"Error: Fastq file {fastq_file_2} is empty for sample {sample}.\n")
    # Print status message
    print(f"Fastq files contain reads for sample {sample}.\n", flush=True)

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
    # Command to run kraken2
    command_kraken2 = [
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
        f"{output_dir}/kraken2/reports/{sample}_kraken2_report.txt"
    ]
    # Run the kraken2 command
    subprocess.run(command_kraken2)
    # Print status message
    print(f"\nFinished kraken2 for sample {sample}\n", flush=True)

def flagstat(sample, output_dir):
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
    # Select the 12th line (amount of properly paired reads)
    twelveth_line = lines[11]
    # Extract the number from the 12th line
    properly_paired = int(twelveth_line.split(' ')[0])
    # Calculate amount of disconcordantly paired reads
    unconcordantly_aligned = total_input - properly_paired
    # Print the results
    print(f"Number of input reads for sample {sample}: {total_input}", flush=True)
    print(f"Number of concordantly aligned reads for sample {sample}: {properly_paired}", flush=True)
    print(f"Number of unconcordantly aligned reads for sample {sample}: {unconcordantly_aligned}", flush=True)
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
    # Open report, loop through lines
    with open(filename, 'r') as file:
        for line in file:
            # Split columns based on tab
            columns = line.split('\t')
            # If homo sapiens is encountered, take value and print
            if 'Homo sapiens' in columns[5]:
                homo_number = int(columns[1])
                print(f"Number of Homo sapiens reads for sample {sample}: {homo_number}", flush=True)
            # If Bacteria is encountered (based on tax id, some species contain 'Bacteria' in name), take value and print
            if re.fullmatch('^2$', columns[4]):
                bac_number = int(columns[1])
                print(f"Number of Bacterial reads for sample {sample}: {bac_number}", flush=True)
    # Warn user if no Homo sapiens classifications
    if homo_number == 0:
        print(f"Warning: No Homo sapiens classifications found in kraken2 report of sample {sample}.\n", flush=True)
    # Warn user if no bacteria classifications found
    if bac_number == 0:
        print(f"Warning: No Bacteria classifications found in kraken2 report of sample {sample}.\n", flush=True)
    # Return values
    return homo_number, bac_number

# Function that creates csv based on outputs of flagstat and extract homo sapiens reads functions
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
    output_file = f"{output_dir}/read_counts.csv"
    # extract homo sapiens reads
    homo_number, bac_number = kraken2_stats(sample, output_dir)
    # Run flagstat for the current sample
    total_input, properly_paired, unconcordantly_aligned = flagstat(sample, output_dir)
    # Write values from each sample to csv
    with open(output_file, 'a', newline='') as outfile:
        out_writer = csv.writer(outfile)
        out_writer.writerow([sample, total_input, properly_paired, unconcordantly_aligned, homo_number, bac_number])

# Function to run kraken-biom
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
    # Run kraken-biom cmd
    subprocess.run(command_kr_biom)
    # Print status message
    print(f"\nCreated {os.path.join(output_dir,'classified.biom')} with taxonomic information of all samples :)\n", flush=True)

# Main function that calls other functions based on samplesheet and passed arguments,
def main():
    """
    Run pipeline with BWA-MEM2 host depletion, samtools, kraken2, biom conversion and read count summarizer.
    Does this by calling run_bwa_mem2, run_samtools, run_kraken2, read_count_extractor functions

    Args:
        --samplesheet (str): Path to CSV file containing sample information.
        --out_dir (str): Path to output directory (not ending with /).

    Returns:
        None
    """
    # Create parser
    parser = argparse.ArgumentParser(description='Run pipeline with BWA-MEM2 host depletion, samtools, kraken2, biom conversion and read count summarizer.')
    # Add the arguments
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to CSV file containing sample information')
    parser.add_argument('--out_dir', type=str, required=True, help='Path to output directory (not ending with /)')
    # Parse!
    args = parser.parse_args()

    # Assign arguments to variables
    samplesheet = args.samplesheet
    output_dir = args.out_dir

    # Create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create subdirectories for bwa-mem2 and kraken2 output if they do not exist
    if not os.path.exists(f"{output_dir}/bwamem2"):
        os.makedirs(f"{output_dir}/bwamem2")
    if not os.path.exists(f"{output_dir}/kraken2"):
        os.makedirs(f"{output_dir}/kraken2")
    if not os.path.exists(f"{output_dir}/kraken2/reports"):
        os.makedirs(f"{output_dir}/kraken2/reports")

    # Create file for summary of extracted readcounts
    output_file = f"{output_dir}/read_counts.csv"

    # Create output file if it does not exist
    with open(output_file, 'w', newline='') as outfile:
        # Create headers
        headers = ['sample', 'total_input_reads', 'human_bwa', 'not_human_bwa', 'human_kraken', 'bacterial_kraken']
        # Write headers to the output file
        out_writer = csv.writer(outfile)
        out_writer.writerow(headers)
    
    # Initialize a counter for the number of samples processed
    num_samples_processed = 0

    # Read sample information from the CSV file and call functions for each sample
    with open(samplesheet, "r") as ssheet:
        ssheet_reader = csv.DictReader(ssheet)
        for row in ssheet_reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            # Run BWA-MEM2 for the current sample
            run_bwa_mem2(sample, fastq_1, fastq_2, output_dir)
            # Run samtools for the current sample
            run_samtools(sample, output_dir)
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
    print(f"\nCreated {output_file} with summarized read counts of all {num_samples_processed} samples :)\n", flush=True)

# If this script is run directly call the main function
if __name__ == '__main__':
    main()
