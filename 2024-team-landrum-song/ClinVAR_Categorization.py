from Bio import Entrez
import os

# Function to search for a gene in ClinVar and save the first record as XML
def search_and_save_gene_info_indiv(gene, email, only_first=True):
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene, retmax=10)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results['IdList']
    print(f"Found IDs for {gene}: {ids}")

    if only_first and ids:
        ids = ids[:1]  # Limit to the first ID
  
    for id in ids:
        # Fetch the ClinVar record for each ID
        fetch_handle = Entrez.esummary(db=db, id=id, rettype="vcv",from_esearch="true")
        xml_data = fetch_handle.read()
        fetch_handle.close()

        # Debug: Print a snippet of the xml_data to ensure it's not empty
        print(f"XML Data snippet for {id}: {xml_data[:500]}")  # Print the first 500 characters

        if not xml_data.strip():  # Check if xml_data is empty or only whitespace
            print(f"No data retrieved for ID {id}.")
            continue  # Skip to the next ID

        # Define the filename and remove if it exists
        filename = f"{gene}_{id}.xml"
        if os.path.exists(filename):
            os.remove(filename)

        # Save the XML data to a file, handling bytes-to-string conversion if necessary
        with open(filename, 'w', encoding='utf-8') as file:
            if isinstance(xml_data, bytes):
                xml_data = xml_data.decode('utf-8')
            file.write(xml_data)
        print(f"Saved {filename}")

# Function to search for a gene in ClinVar and save the first record as XML
def search_and_save_gene_info_combined(gene, email, only_first=True):
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene, retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results['IdList']
    print(f"Found IDs for {gene}: {ids}")

    if only_first and ids:
        ids = ids[:1]  # Limit to the first ID
  
  
    # Fetch the ClinVar record for each ID
    id_list = ",".join(ids)
    fetch_handle = Entrez.efetch(db=db, id=id_list, rettype="vcv",from_esearch="true")
    xml_data = fetch_handle.read()
    fetch_handle.close()
    
    # Define the filename and remove if it exists
    filename = f"{gene}.xml"
    if os.path.exists(filename):
        os.remove(filename)

    # Save the XML data to a file, handling bytes-to-string conversion if necessary
    with open(filename, 'w', encoding='utf-8') as file:
        if isinstance(xml_data, bytes):
            xml_data = xml_data.decode('utf-8')
        file.write(xml_data)
    print(f"Saved {filename}")

# Example usage
# genes = ["LDLR", "KCNQ1", "251590", "USH2A", "SCN5A", "TSC1"]
genes = ["LDLR"]
email = "wengang.zhang@nih.gov"

# To only search one gene and the first ID, set only_first to True
only_first = False

for gene in genes:
    search_and_save_gene_info_combined(gene, email, only_first=only_first)
