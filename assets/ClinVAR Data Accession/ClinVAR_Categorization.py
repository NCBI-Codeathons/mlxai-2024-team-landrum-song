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
    search_handle = Entrez.esearch(db=db, term=gene + "[Gene Name]", retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]  # ['250925']#
    # print(f"Found IDs for {gene}: {ids}")

    print(
        f"Number of IDs found for {gene}: {len(ids)}"
    )  # print the number of IDs found

    if only_first and ids:
        ids = ids[:1]  # Limit to the first ID

    # Fetch the ClinVar record for each ID
    id_list = ",".join(ids)
    fetch_handle = Entrez.esummary(
        db=db, id=id_list, rettype="vcv", from_esearch="true"
    )
    xml_data = fetch_handle.read()
    fetch_handle.close()

    # Load and parse the XML file
    root = ET.fromstring(xml_data)

    # Initialize a list to store extracted data
    data = []

    # Extract relevant information
    for doc_summary in root.findall(".//DocumentSummary"):
        variant_id = doc_summary.get("uid")
        rcv_ids = [rcv.text for rcv in doc_summary.findall(".//rcv/string")]
        medgen_entries = [
            (xref.find(".//db_source").text, xref.find(".//db_id").text)
            for xref in doc_summary.findall(".//trait_xref")
            if xref.find(".//db_source").text == "MedGen"
        ]
        trait_names = [name.text for name in doc_summary.findall(".//trait_name")]

        # Assuming each RCV ID corresponds to a trait name and MedGen entry in order
        for i, rcv_id in enumerate(rcv_ids):
            medgen_db, medgen_id = (
                medgen_entries[i] if i < len(medgen_entries) else ("Unknown", "Unknown")
            )
            trait_name = trait_names[i] if i < len(trait_names) else "Unknown"

            # Construct a dictionary for each variant's information
            variant_info = {
                "variant_id": variant_id,
                "RCV_id": rcv_id,
                "MedGen": {"db": medgen_db, "id": medgen_id},
                "trait_name": trait_name,
            }

            data.append(variant_info)

    # Display the structured data extracted
    data[:10]  # Display the first 10 entries for brevity

    # Save to XMLof the variant id, condition
    xml_filename = gene + "_variants_extracted.xml"
    # Create the root element
    root = Element("Variants")

    # Iterate through each item in the data list and create XML structure
    for item in data:
        variant_el = SubElement(root, "Variant", id=item["variant_id"])
        RCV_el = SubElement(variant_el, "RCV", id=item["RCV_id"])
        # MedGen_el = SubElement(
        #     RCV_el, "MedGen", db=item["MedGen"]["db"], id=item["MedGen"]["id"]
        # )
        trait_name_el = SubElement(RCV_el, "TraitName")
        trait_name_el.text = item["trait_name"]

    # Create an ElementTree object from the root element
    tree = ElementTree(root)
    if os.path.exists(xml_filename):
        os.remove(xml_filename)
    tree.write(xml_filename, encoding="UTF-8", xml_declaration=True)
    print(f"Saved variant information to {xml_filename}")

    # Define the filename and remove if it exists
    filename = f"{gene}.xml"
    if os.path.exists(filename):
        os.remove(filename)

    # Save the XML data to a file, handling bytes-to-string conversion if necessary
    with open(filename, "w", encoding="utf-8") as file:
        if isinstance(xml_data, bytes):
            xml_data = xml_data.decode("utf-8")
        file.write(xml_data)
    print(f"Saved {filename}")

    # print("Visualizing the first 10 entries in JSON format:")
    # print(json.dumps(data[:10], indent=4))

    # Save the data to a JSON file
    json_filename = gene + "_variants_extracted.json"
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)
    print(f"Saved variant information to {json_filename}")


# Example usage
genes = ["LDLR"]
email = "bkkesler4@gmail.com"

# To only search one gene and the first ID, set only_first to True
only_first = False

for gene in genes:
    search_and_save_gene_info_combined(gene, email, only_first=only_first)
