# Script to run samtools flagstat on all SAM files in a directory and saves all outputs to tsv file
# Make sure samtools is installed in your environment
import os
import subprocess
import argparse

# Create parser
parser = argparse.ArgumentParser(
    description=('Run samtools flagstat on all SAM files in a directory and save all outputs to a TSV file.')
)

# Add arguments
parser.add_argument('--sam_dir', '-d', type=str, required=True,
                    help='Path to directory containing SAM files')
parser.add_argument('--output_file', '-o', type=str, required=True,
                    help='Path to output TSV file')

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
