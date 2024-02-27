#!/usr/bin/env python3

"""
This script outputs diseases related to genetic variations of LDLR
"""

import os
import xml.etree.ElementTree as ET

from Bio import Entrez


def search_and_save_gene_info_combined(gene: str, email: str, only_first: bool) -> None:
    """
    TODO
    """
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene, retmax=10000)
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

    # Extract ClassifiedCondition elements
    classified_conditions = [
        condition.text
        for condition in root.findall(".//ClassifiedCondition")
        if condition.text
    ]

    # Define the filename for the classified conditions
    classified_filename = f"{gene}_classified_conditions.txt"
    if os.path.exists(classified_filename):
        os.remove(classified_filename)

    # Save the ClassifiedCondition data to a file
    with open(classified_filename, "w", encoding="utf-8") as file:
        for condition in classified_conditions:
            file.write(condition + "\n")
    print(f"Saved classified conditions to {classified_filename}")


def main() -> None:
    """
    TODO
    """

    # Example usage parameters
    genes = ["LDLR", "251590"]  # ["LDLR", "KCNQ1", "USH2A", "SCN5A", "TSC1", "251590"]
    email = "wengang.zhang@nih.gov"
    only_first = False  # To only search one gene and the first ID for testing set only_first to True

    for gene in genes:
        search_and_save_gene_info_combined(gene, email, only_first=only_first)


if __name__ == "__main__":
    main()
