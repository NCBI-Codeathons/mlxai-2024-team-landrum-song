from Bio import Entrez
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, ElementTree
import os
import json


# Function to search for a gene in ClinVar and save the first record as XML
def search_and_save_gene_info_combined(gene, email, only_first=True):
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene+"[Gene Name]"+"AND \"single gene\"[prop]", retmax=10000) #
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results['IdList'] #['250925']#
    #print(f"Found IDs for {gene}: {ids}")
    #print(f"Number of IDs found for {gene}: {len(ids)}") #print the number of IDs found

    # Fetch the ClinVar record for each ID
    id_list = ",".join(ids)

    total_records = len(ids)
    print(f"Total records found for {gene}: {total_records}")

    retmax = 1000  # Number of records to fetch per request
    fetched_records = 0

    data = []  # Initialize a list to store fetched data

    while fetched_records < total_records:
        fetch_handle = Entrez.efetch(db=db, id=id_list, rettype="vcv",  retstart=fetched_records, retmax=retmax, from_esearch="true") #       
        try:
            xml_data = fetch_handle.read()
            print("fetch complete")
        except IncompleteRead as e:
            xml_data = e.partial  # Use partial data if an IncompleteRead exception occurs
            print("fetch incomplete")
        finally:
            fetch_handle.close()

        # Process the XML data
        root = ET.fromstring(xml_data)
        # Extract and process the relevant information from `root` as per your existing logic
        # Extract relevant information
        for variant in root.findall(".//VariationArchive"):
            variant_id = variant.get("VariationID") or "Unknown"
            for rcv_accession in variant.findall(".//RCVAccession"):
                rcv_id = rcv_accession.get("Accession") or "Unknown"
                for classified_condition in rcv_accession.findall(".//ClassifiedCondition"):
                    condition_text = classified_condition.text or "Unknown"
                    condition_db = classified_condition.get("DB") or "Unknown"
                    condition_id = classified_condition.get("ID") or "Unknown"
                    data.append({
                        "variant_id": variant_id,
                        "RCV": rcv_id,
                        "condition": {
                            "db": condition_db,
                            "id": condition_id,
                            "text": condition_text,
                        }
                    })
        # Update fetched_records to reflect the number of records processed
        #print(data[:5])  # Display the first 10 entries for brevity
        fetched_records += retmax
        print(fetched_records)
 
    # Display the structured data extracted
    print(data[:10])  # Display the first 10 entries for brevity

    # Save to XMLof the variant id, condition
    xml_filename = gene + "_variants_extracted.xml"
  
    # Save to XML of the variant id, condition
    xml_filename = gene + "_variants_extracted.xml"
    root_el = ET.Element("Variants")
    for variant in data:
        variant_el = ET.SubElement(root_el, "Variant")
        ET.SubElement(variant_el, "VariantID").text = variant["variant_id"]
        ET.SubElement(variant_el, "RCV").text = variant["RCV"]
        condition_el = ET.SubElement(variant_el, "ClassifiedCondition")
        condition_el.set("DB", variant["condition"]["db"])
        condition_el.set("ID", variant["condition"]["id"])
        condition_el.text = variant["condition"]["text"]
    tree = ET.ElementTree(root_el)

    # Create an ElementTree object from the root element
    tree = ElementTree(root)
    if os.path.exists(xml_filename):
        os.remove(xml_filename)
    tree.write(xml_filename, encoding='UTF-8', xml_declaration=True)
    print(f"Saved variant information to {xml_filename}")


    
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


    #print("Visualizing the first 10 entries in JSON format:")
    #print(json.dumps(data[:10], indent=4))

    # Save the data to a JSON file
    json_filename = gene + "_variants_extracted.json"
    with open(json_filename, 'w', encoding='utf-8') as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)
    print(f"Saved variant information to {json_filename}")

    # Extract unique disease names
    unique_disease_names = {item['condition']['text'] for item in data}

    # Save unique disease names to JSON
    unique_disease_json_filename = gene + "_unique_diseases.json"

    with open(unique_disease_json_filename, 'w', encoding='utf-8') as json_file:
        # Convert the set to a list for JSON serialization
        json.dump(list(unique_disease_names), json_file, ensure_ascii=False, indent=4)
    print(f"Saved unique disease names to {unique_disease_json_filename}")

    # Assigning an ID to each unique condition and formatting the data
    formatted_conditions = [{"cid": idx + 1, "condition": condition_text} for idx, condition_text in enumerate(unique_disease_names)]

    # Saving the formatted unique condition texts to a JSON file
    json_filename = gene + "_formatted_unique_conditions.json"
    with open(json_filename, 'w', encoding='utf-8') as json_file:
        json.dump(formatted_conditions, json_file, ensure_ascii=False, indent=4)

    print(f"Saved formatted unique condition texts to {json_filename}")


# Example usage
genes = ["LDLR","USH2A", "TSC1", "KCNQ1", "SCN5A"]
email = "wengang.zh@gmail.com"

# To only search one gene and the first ID, set only_first to True
only_first = False

for gene in genes:
    search_and_save_gene_info_combined(gene, email, only_first=only_first)
