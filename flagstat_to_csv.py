# Script to run samtools flagstat on all SAM files in a directory and saves all outputs to tsv file
# Also parses output to separate file with only filename and total reads
# Make sure samtools is installed in your environment
import os
import subprocess
import argparse
import csv

# Create parser
parser = argparse.ArgumentParser(
    description=('Run samtools flagstat on all BAM files in a directory and save all outputs to a TSV file.')
)

# Add arguments
parser.add_argument('--sam_dir', '-d', type=str, required=True,
                    help='Path to directory containing BAM files')
parser.add_argument('--output_file', '-o', type=str, required=True,
                    help='Path to output file')

# Parse! :)
args = parser.parse_args()

# Get list of all SAM files in the directory
sam_files = [file for file in os.listdir(args.sam_dir) if file.endswith('.bam')]

# Run samtools flagstat on the remaining SAM files and append the output (excluding header) to the TSV file
for sam_file in sam_files:
    sam_path = os.path.join(args.sam_dir, sam_file)
    # Print sam_file in tsv file
    with open(args.output_file, 'a') as f:
        f.write(f'{sam_file}\n')
    subprocess.run(f'samtools flagstat -@ 32 -O tsv {sam_path} >> {args.output_file}', shell=True)

# Initialize lists to store filenames and total reads
filenames = []
total_reads = []

# Open and read the TSV file line by line
with open(args.output_file, "r") as file:
    reader = csv.reader(file, delimiter="\t")
    for row in reader:
        # Check if the line has a single element ending with '.bam' (indicating it's a filename)
        if len(row) == 1 and row[0].endswith(".bam"):
            filename = row[0]
            filenames.append(filename)
        # Check if the line has "total (QC-passed reads + QC-failed reads)" to get total reads
        elif len(row) > 2 and "total (QC-passed reads + QC-failed reads)" in row[2]:
            total_reads.append(row[0])

# Combine filenames and total reads into a list of tuples
data = list(zip(filenames, total_reads))

# Create filename based on input but parsed
parsed_output_file = args.output_file.split('.')[0] + '_parsed.csv'

# Save the result to a new CSV file
with open(parsed_output_file, "w", newline="") as output_file:
    writer = csv.writer(output_file)
    writer.writerow(["Filename", "Total Reads"])  # Write header
    writer.writerows(data)  # Write data rows

# # Print the result for confirmation
# for filename, total in data:
#     print(f"{filename} | {total}")
