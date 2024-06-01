# Script that runs Bracken default on strain, species and genus levels on all kreport files in a dir
# Make sure Bracken is installed in your environment and the Bracken database is build
import os
import subprocess
import argparse

# Create parser
parser = argparse.ArgumentParser(
    description=('Run Bracken on all Kraken2 reports in a directory.')
)
# Add arguments
parser.add_argument('--kr_dir', '-k', type=str, required=True,
                    help='Path to directory containing Kraken2 reports')
parser.add_argument('--db_path', '-d', type=str, required=True,
                    help='Path to Bracken database')
parser.add_argument('--output_dir', '-o', type=str, required=True,
                    help = 'Path to output directory')
# Parse
args = parser.parse_args()

# Output and report directories
output_dir_reports = os.path.join(args.output_dir, 'bracken_reports')

# Make dirs if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(output_dir_reports):
    os.makedirs(output_dir_reports)

# Iterate over files
for filename in os.listdir(args.kr_dir):
    if filename.endswith("kraken2_report.txt"):
        # Full krakenrrport path
        kreport = os.path.join(args.kr_dir, filename)
        
        # Full bracken kreport path species
        br_kreport_S = os.path.join(output_dir_reports, filename.replace("kraken2_report.txt", "bracken_report_S.txt"))
        # Full bracken default output path species
        br_output_S = os.path.join(args.output_dir, filename.replace("kraken2_report.txt", "bracken_out_S.txt"))

        # Full bracken kreport path genus
        br_kreport_G = os.path.join(output_dir_reports, filename.replace("kraken2_report.txt", "bracken_report_G.txt"))
        # Full bracken default output path genus
        br_output_G = os.path.join(args.output_dir, filename.replace("kraken2_report.txt", "bracken_out_G.txt"))

        # Full bracken kreport path strain
        br_kreport_S1 = os.path.join(output_dir_reports, filename.replace("kraken2_report.txt", "bracken_report_S1.txt"))
        # Full bracken default output path genus
        br_output_S1 = os.path.join(args.output_dir, filename.replace("kraken2_report.txt", "bracken_out_S1.txt"))
        
        # Run bracken genus
        subprocess.run(["bracken", "-d", args.db_path, "-i", kreport, "-o", br_output_G, "-w", br_kreport_G, "-l", "G"])

        # Run bracken species
        subprocess.run(["bracken", "-d", args.db_path, "-i", kreport, "-o", br_output_S, "-w", br_kreport_S, "-l", "S"])

        # Run bracken strain
        subprocess.run(["bracken", "-d", args.db_path, "-i", kreport, "-o", br_output_S1, "-w", br_kreport_S1, "-l", "S1"])