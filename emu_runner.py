import os
import subprocess
import argparse
import csv

# Create parser
parser = argparse.ArgumentParser(
    description=('Run emu abundance command based on a samplesheet.')
)
# Add arguments
parser.add_argument('--samplesheet', '-s', type=str, required=True,
                    help='Path to the samplesheet CSV file')
parser.add_argument('--output_dir', '-o', type=str, required=True,
                    help='Path to directory where output will be saved')
parser.add_argument('--db', '-d', type=str, default="/mnt/SCS/data_2/Users/Birgit/db_index/emu/",
                    help='Path to the database index for emu')
parser.add_argument('--threads', '-t', type=int, default=16,
                    help='Number of threads to use (default: 16)')
parser.add_argument('--type', '-y', type=str, default='map-ont',
                    help='Type of sequencing data (default: map-ont)')

# Parse arguments
args = parser.parse_args()

# Read the samplesheet (2 columns, sample (name) and fastq_1 (path))
with open(args.samplesheet, 'r') as ssheet:
    reader = csv.DictReader(ssheet)
    # Loop over rows in the samplesheet
    for row in reader:
        input_path = row['fastq_1']
        # Create output path based on sample name (first column)
        output_path = os.path.join(args.output_dir, row['sample'])
        # Create emu command
        command = [
            'emu', 'abundance', input_path,
            '--type', args.type,
            '--db', args.db,
            '--output-dir', output_path,
            '--threads', str(args.threads)
        ]
        # Run the command
        subprocess.run(command)