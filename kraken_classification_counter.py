# This script counts classifications in a Kraken2 output file (not `--report`) with scientific names (`--use-names`)

# Import modules
import csv
import argparse

# Function to check kraken2 output format, exit if format is not as expected
def filechecker(kr2_output):
    with open(kr2_output, 'r') as file:
        first_line = file.readline()
        parts = first_line.split("\t")
        # Check if file contains 5 columns
        if len(parts) != 5:
            print("Error: File format is not as expected.")
            return False
        # Check if first column contains C/U (classified/unclassified)
        if not parts[0] in ['C', 'U']:
            print("Error: First column should be 'C' or 'U'.")
            return False
        # Check if third column contains "taxid"
        if not "taxid" in parts[2]:
            print("Error: Third column should contain 'taxid'.")
            return False
        # Check if fourth column contains "|"
        if not "|" in parts[3]:
            print("Error: Fourth column should contain '|'.")
            return False
        # Check if fifth column contains ":"
        if not ":" in parts[4]:
            print("Error: Fifth column should contain ':'.")
            return False
    return True


def counter(kr2_output, extra_classification, human_only, output_other):
    # Intitialize counts
    read_count = 0
    unclass_count = 0
    root_count = 0
    homo_sapiens_count = 0
    extra_count = 0
    # Intialize dictionary
    other_classifications = {}
    # Set for other classifications
    other_classifications = set()

    # Open tabular file to read
    with open(kr2_output, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        # Iterate through rows of file
        if human_only:
            for row in reader:
                # Increase counters when specific classification is encountered
                read_count += 1
                if row[2] == 'Homo sapiens (taxid 9606)':
                    homo_sapiens_count += 1
        elif output_other:
            for row in reader:
                # Increase counters when specific classification is encountered
                read_count += 1
                if row[2] == 'Homo sapiens (taxid 9606)':
                    homo_sapiens_count += 1
                elif output_other and row[2] == extra_classification:
                    extra_count += 1
                elif row[2] == 'unclassified (taxid 0)':
                    unclass_count += 1
                elif row[2] == 'root (taxid 1)':
                    root_count += 1
                # If not human, unclassified, root or extra, save classification in dictonary with occurence count
                else:
                    if row[2] in other_classifications:
                        other_classifications[row[2]] +=1
                    else:
                        other_classifications[row[2]] = 1
        else:
            # Increase counters when specific classification is encountered
            for row in reader:
                read_count += 1
                if row[2] == 'Homo sapiens (taxid 9606)':
                    homo_sapiens_count += 1
                elif extra_classification and row[2] == extra_classification:
                    extra_count += 1
                elif row[2] == 'unclassified (taxid 0)':
                    unclass_count += 1
                elif row[2] == 'root (taxid 1)':
                    root_count += 1
    return read_count, homo_sapiens_count, extra_count, unclass_count, root_count

# Function that prints the found counts by counter function
def printer(kr2_output, extra_classification, human_only, output_other):
    # Call the counter function and get the counts
    read_count, homo_sapiens_count, extra_count, unclass_count, root_count = counter()
    # Only return human counts if specified by user
    if human_only:
        print(f"Homo sapiens classifications: {homo_sapiens_count}")
        if output_other:
            print(f"Error: Since '--human-only' was specified, other classifications were not saved.")

    # Return all counts
    else:
        print(f"Total input reads: {read_count}")
        print(f"Homo sapiens classifications: {homo_sapiens_count}")
        print(f"Root classifications: {root_count}")
        print(f"Unclassified: {unclass_count}")
        # Return extra classification count if specified
        if extra_classification:
            print(f"{extra_classification} classifications: {extra_count}")

    # Return other classifications if specified by user
    if output_other:
        if human_only:
            print(f"Error: Since '--output_other' was specified, '--human_only' does not work.")
        else:
            print(f"Other classifications:")
            for classification, count in other_classifications.items():
                print(f"{classification}: {count}")

# Main function with argparse
def main():
    # Create parser
    parser = argparse.ArgumentParser(description='Count classifications in a Kraken2 report file with scientific names.')

    # Create argparse thingy
    parser = argparse.ArgumentParser(
        prog='Kraken2 classification summarizer',
        description='Count classifications in a Kraken2 report file with scientific names.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='Thanks for using this tool! :)'
    )

    # Add arguments
    parser.add_argument('kr_out_path', type=str,
                        help = 'Path to the kraken2 output file.')
    parser.add_argument('-e', '--extra', dest='extra_classification', type=str, default=None,
                        help = 'Extra classification that you want to count.')
    parser.add_argument('--human_only', action='store_true',
                        help = 'Only print the amount of Homo sapiens classifications.')
    parser.add_argument('--output_other', action='store_true',
                        help = 'Show a list of all unique classifications that are not human, root or unclassified.' )

    # Parse! :)
    args = parser.parse_args()

    # Call functions
    filechecker(args.kr_out_path)
    #counter(args.kr_out_path, args.extra, args.human_only, args.output_other)
    #Printer? Or integrate into counter

if __name__ == "__main__":
    main()