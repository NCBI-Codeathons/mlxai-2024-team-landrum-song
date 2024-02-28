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
import http
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple
from xml.etree.ElementTree import ElementTree

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


async def pull_and_extract_clinvar_data(gene: str, email: str) -> List[Dict]:
    """
    Download the complete ClinVar XML dataset for the gene of interest.
    """
    Entrez.email = email
    db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=db, term=gene + "[Gene Name]", retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]

    # Fetch the ClinVar record for each ID
    id_list = ",".join(ids)

    total_records = len(ids)
    rprint(f"Total records found for {gene}: {total_records}")

    retmax = 1000  # Number of records to fetch per request
    fetched_records = 0

    data = []  # Initialize a list to store fetched data

    while fetched_records < total_records:
        fetch_handle = Entrez.efetch(
            db=db,
            id=id_list,
            rettype="vcv",
            retstart=fetched_records,
            retmax=retmax,
            from_esearch="true",
        )
        try:
            xml_data = fetch_handle.read()
        except http.client.IncompleteRead as e:
            xml_data = (
                e.partial
            )  # Use partial data if an IncompleteRead exception occurs
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
                for classified_condition in rcv_accession.findall(
                    ".//ClassifiedCondition"
                ):
                    condition_text = classified_condition.text or "Unknown"
                    condition_db = classified_condition.get("DB") or "Unknown"
                    condition_id = classified_condition.get("ID") or "Unknown"
                    data.append(
                        {
                            "variant_id": variant_id,
                            "RCV": rcv_id,
                            "condition": {
                                "db": condition_db,
                                "id": condition_id,
                                "text": condition_text,
                            },
                        }
                    )
        # Update fetched_records to reflect the number of records processed
        fetched_records += retmax

    condition_set = set(record.get("condition").get("text") for record in data)
    rprint(f"{len(condition_set)} unique conditions found for gene {gene}.")

    return data


async def save_variant_xml(gene: str, data: List) -> None:
    """
    Save the extracted variant XML data for the purposes of this codeathon.
    """

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

    # Create an ElementTree object from the root element
    tree = ElementTree(root_el)
    if os.path.exists(xml_filename):
        os.remove(xml_filename)
    tree.write(xml_filename, encoding="UTF-8", xml_declaration=True)


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

    # collect the list of genes
    with open(gene_file, "r", encoding="utf8") as input_handle:
        genes = [gene.strip() for gene in input_handle]

    # process each gene asynchronously
    for gene in genes:
        extracted_data = await pull_and_extract_clinvar_data(gene, email)
        await save_variant_xml(gene, extracted_data)
        await save_variant_json(gene, extracted_data)
        rprint(f"Data retrieval, extraction, and writing complete for {gene}")


if __name__ == "__main__":
    asyncio.run(main())
