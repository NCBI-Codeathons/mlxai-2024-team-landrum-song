from Bio import Entrez
import os
import xml.etree.ElementTree as ET


def search_and_save_gene_info_combined(gene, email, only_first=True):
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene, retmax=100)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]
    print(f"Found IDs for {gene}: {ids}")

    if only_first and ids:
        ids = ids[:1]  # Limit to the first ID

    # Fetch the ClinVar records for the IDs
    id_list = ",".join(ids)
    fetch_handle = Entrez.efetch(db=db, id=id_list, rettype="vcv", retmode="xml")
    xml_data = fetch_handle.read()
    fetch_handle.close()

    # Parse the XML data
    root = ET.fromstring(xml_data)

    # Initialize a list to store ClassifiedCondition data
    classified_conditions = []

    # Extract ClassifiedCondition elements
    for condition in root.findall(".//ClassifiedCondition"):
        condition_text = condition.text
        if condition_text:
            classified_conditions.append(condition_text)

    # Define the filename for the classified conditions
    classified_filename = f"{gene}_classified_conditions.txt"
    if os.path.exists(classified_filename):
        os.remove(classified_filename)

    # Save the ClassifiedCondition data to a file
    with open(classified_filename, "w", encoding="utf-8") as file:
        for condition in classified_conditions:
            file.write(condition + "\n")
    print(f"Saved classified conditions to {classified_filename}")


# Example usage parameters
genes = ["LDLR", "KCNQ1", "USH2A", "SCN5A", "TSC1"]
email = "bkkesler4@gmail.com"

# To only search one gene and the first ID, set only_first to True
only_first = False

for gene in genes:
    search_and_save_gene_info_combined(gene, email, only_first=only_first)
