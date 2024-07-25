# Script that extracts read counts from all fastp output json files and puts results in csv format
import os
import json
import csv

# Directory with fastp outputs
fastp_out_dir = "/mnt/FS2/data_2/Users/Birgit/nf_out/CRC_60_samples/fastp/"

# Get a list of all json files in the directory
json_files = [f for f in os.listdir(fastp_out_dir) if f.endswith('.json')]

# Initialize two empty lists
before_filtering = []
after_filtering = []

# Loop through each file in the list of json files
for file in json_files:
    # Open each file and load the json data
    with open(os.path.join(fastp_out_dir, file)) as f:
        data = json.load(f)
    
    # Extract the "total_reads" value from the "before_filtering" and "after_filtering" sections
    before_filtering.append(data['summary']['before_filtering']['total_reads'])
    after_filtering.append(data['summary']['after_filtering']['total_reads'])

# Create a csv file and write the sample name, before_filtering and after_filtering values
with open('output.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['sample', 'before_filtering', 'after_filtering'])
    
    for file, before, after in zip(json_files, before_filtering, after_filtering):
        writer.writerow([os.path.splitext(file)[0], before, after])