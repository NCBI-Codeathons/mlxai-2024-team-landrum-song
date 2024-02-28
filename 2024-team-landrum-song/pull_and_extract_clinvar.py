#!/usr/bin/env python3

"""
This script interfaces with NCBI ClinVar via the Entrez functionality included
in Biopython to download XML-formatted variant data, extract information about
the condition(s) associated with each variant, and optionally save that
information as a JSON and/or a simplified XML.

```
usage: pull_and_extract_clinvar.py [-h] --genelist GENELIST --email EMAIL

options:
  -h, --help            show this help message and exit
  --genelist GENELIST, -g GENELIST
                        Text file with one, headerless column of one gene name per line.
  --email EMAIL, -e EMAIL
                        Email for NCBI ClinVar tracking purposes.
```
"""

import argparse
import asyncio
import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple
from xml.etree.ElementTree import Element, ElementTree, SubElement

from Bio import Entrez
from rich import print as rprint


def parse_command_line_args() -> Tuple[Path, str]:
    """
    Parse command line arguments, returning a file that lists genes of
    interest and an email for NCBI tracking purposes.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genelist",
        "-g",
        type=Path,
        required=True,
        help="Text file with one, headerless column of one gene name per line.",
    )
    parser.add_argument(
        "--email",
        "-e",
        type=str,
        required=True,
        help="Email for NCBI ClinVar tracking purposes.",
    )
    args = parser.parse_args()

    return args.genelist, args.email


async def pull_clinvar_xml(gene: str, email: str, only_first: bool) -> str:
    """
    Download the complete ClinVar XML dataset for the gene of interest.
    """
    Entrez.email = email
    ncbi_db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=ncbi_db, term=gene, retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]

    rprint(f"Number of IDs found for {gene}: {len(ids)}")

    if only_first and ids:
        ids = ids[:1]  # Limit to the first ID

    # Fetch the ClinVar record for each ID
    id_list = ",".join(ids)
    fetch_handle = Entrez.esummary(
        db=ncbi_db, id=id_list, rettype="vcv", from_esearch="true"
    )
    xml_data = fetch_handle.read()
    fetch_handle.close()

    # write out raw xml_data
    await save_full_xml(gene, xml_data)

    # Load and parse the XML file
    root = ET.fromstring(xml_data)

    return root


async def extract_variant_data(root) -> List[Dict]:
    """
    Extract RCV IDs, MedGen IDs, and other metadata for each variant and
    return that information as a list of dictionaries.
    """

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

    return data


async def save_variant_xml(gene: str, data: List) -> None:
    """
    Save the extracted variant XML data for the purposes of this codeathon.
    """
    # Save to XMLof the variant id, condition
    xml_filename = gene + "_variants_extracted.xml"
    # Create the root element
    new_root = Element("Variants")

    # Iterate through each item in the data list and create XML structure
    for item in data:
        variant_el = SubElement(new_root, "Variant", id=item["variant_id"])
        rcv_el = SubElement(variant_el, "RCV", id=item["RCV_id"])
        _medgen_el = SubElement(
            rcv_el, "MedGen", db=item["MedGen"]["db"], id=item["MedGen"]["id"]
        )
        trait_name_el = SubElement(rcv_el, "TraitName")
        trait_name_el.text = item["trait_name"]

    # Create an ElementTree object from the root element
    tree = ElementTree(new_root)
    if os.path.exists(xml_filename):
        os.remove(xml_filename)
    tree.write(xml_filename, encoding="UTF-8", xml_declaration=True)
    rprint(f"Saved variant information to {xml_filename}")


async def save_full_xml(gene: str, xml_data) -> None:
    """
    Save the raw XML data as an XML file for future usage for other
    purposes.
    """

    # Define the filename and remove if it exists
    filename = f"{gene}.xml"
    if os.path.exists(filename):
        os.remove(filename)

    # Save the XML data to a file, handling bytes-to-string conversion if necessary
    with open(filename, "w", encoding="utf-8") as file:
        if isinstance(xml_data, bytes):
            new_xml_data = xml_data.decode("utf-8")
            file.write(new_xml_data)
        else:
            file.write(xml_data)
    rprint(f"Saved {filename}")


async def save_variant_json(gene: str, data) -> None:
    """
    Save the extracted variant data in JSON format for the purposes of this codeathon.
    """
    # Save the data to a JSON file
    json_filename = gene + "_variants_extracted.json"
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)
    rprint(f"Saved variant information to {json_filename}")


async def main() -> None:
    """
    Coordinate the flow of data through the above functions within an
    asynchronous runtime that separates network/IO-bound tasks from
    CPU-bound tasks.
    """

    # Example usage parameters
    gene_file, email = parse_command_line_args()

    # To only search one gene and the first ID for testing set only_first to True
    only_first = False

    # collect the list of genes
    with open(gene_file, "r", encoding="utf8") as input_handle:
        genes = [gene.strip() for gene in input_handle]

    # process each gene asynchronously
    for gene in genes:
        xml_data = await pull_clinvar_xml(gene, email, only_first=only_first)
        extracted_data = await extract_variant_data(xml_data)
        await save_variant_xml(gene, extracted_data)
        await save_variant_json(gene, extracted_data)
        rprint(f"Data retrieval, extraction, and writing complete for {gene}")


if __name__ == "__main__":
    asyncio.run(main())
