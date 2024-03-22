# This script fetches the NCBI taxonomy database IDs for a list of bacteria names, and puts them in a csv file :)

import requests
import xml.etree.ElementTree as ET

# List of bacteria names (from paper Salter et al.)
    # Modify if needed
data = "Afipia,Aquabacterium,Asticcacaulis,Aurantimonas,Beijerinckia,Bosea,Bradyrhizobium,Brevundimonas,Caulobacter,Craurococcus,Devosia,Hoeflea,Mesorhizobium,Methylobacterium,Novosphingobium,Ochrobactrum,Paracoccus,Pedomicrobium,Phyllobacterium,Rhizobium,Roseomonas,Sphingobium,Sphingomonas,Sphingopyxis,Acidovorax,Azoarcus,Azospira,Burkholderia,Comamonas,Cupriavidus,Curvibacter,Delftia,Duganella,Herbaspirillum,Janthinobacterium,Kingella,Leptothrix,Limnobacter,Massilia,Methylophilus,Methyloversatilis,Oxalobacter,Pelomonas,Polaromonas,Ralstonia,Schlegelella,Sulfuritalea,Undibacterium,Variovorax,Acinetobacter,Enhydrobacter,Enterobacter,Escherichia,Nevskia,Pseudomonas,Pseudoxanthomonas,Psychrobacter,Stenotrophomonas,Xanthomonas,Aeromicrobium,Arthrobacter,Beutenbergia,Brevibacterium,Corynebacterium,Curtobacterium,Dietzia,Geodermatophilus,Janibacter,Kocuria,Microbacterium,Micrococcus,Microlunatus,Patulibacter,Propionibacterium,Rhodococcus,Tsukamurella,Abiotrophia,Bacillus,Brevibacillus,Brochothrix,Facklamia,Paenibacillus,Streptococcus,Chryseobacterium,Dyadobacter,Flavobacterium,Hydrotalea,Niastella,Olivibacter,Pedobacter,Wautersiella,Deinococcus"

# Path to output csv file
    # Modify if needed
output_path = "bacteria_ids.csv"

# Base URL for the NCBI taxonomy database search with eutils
base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term="

# Split the data into a list of bacteria names
contaminants = data.split(",")

# Create a dictionary to store the fetched ids
ids = {}

# Loop through the contaminant list
for bac in contaminants:
    # Send request to NCBI taxonomy database based on base url + name
    response = requests.get(base_url + bac)
    
    # Parse the XML response
    root = ET.fromstring(response.content)
    
    # Extract the ID and store it in the dictionary
    ids[bac] = root.find('IdList').find('Id').text

# Write dictonary to csv file
with open(output_path, 'w') as f:
    # Header
    f.write("Bacteria,NCBI ID\n")
    # Write each combination to new row in csv file
    for key in ids.keys():
        f.write("%s,%s\n"%(key,ids[key]))