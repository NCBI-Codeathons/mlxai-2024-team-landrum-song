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
import http
import json
import os
import urllib
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from xml.etree.ElementTree import ElementTree

from Bio import Entrez
from rich import print as rprint


def parse_command_line_args() -> Tuple[Path, Optional[str], str]:
    """
    Parse command line arguments, returning a file that lists genes of
    interest and an email for NCBI tracking purposes.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genelist",
        "-G",
        type=Path,
        required=False,
        help="Text file with one, headerless column of one gene name per line.",
    )
    parser.add_argument(
        "--gene",
        "-g",
        type=str,
        required=False,
        default=None,
        help="Single gene to query ClinVar for.",
    )
    parser.add_argument(
        "--email",
        "-e",
        type=str,
        required=True,
        help="Email for NCBI ClinVar tracking purposes.",
    )
    args = parser.parse_args()

    return args.genelist, args.gene, args.email


async def pull_clinvar_data(
    db: str, id_list: List[str], fetched_records: int, retmax: int
):
    """
    TODO
    """
    return Entrez.efetch(
        db=db,
        id=id_list,
        rettype="vcv",
        retstart=fetched_records,
        retmax=retmax,
        from_esearch="true",
    )


async def pull_and_extract_data(gene: str, email: str) -> List[Dict]:
    """
    Download the complete ClinVar XML dataset for the gene of interest.
    """
    Entrez.email = email
    ncbi_db = "clinvar"

    # Search for the gene in ClinVar
    search_handle = Entrez.esearch(db=ncbi_db, term=gene + "[Gene Name]", retmax=10000)
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
        fetch_handle = await pull_clinvar_data(
            ncbi_db, id_list, fetched_records, retmax
        )
        try:
            xml_data = fetch_handle.read()
        except http.client.IncompleteRead as incomplete:
            xml_data = (
                incomplete.partial
            )  # Use partial data if an IncompleteRead exception occurs
        except urllib.error.HTTPError as _:
            continue
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


async def save_variant_xml(gene: str, data: List) -> str:
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

    return xml_filename


async def save_variant_json(gene: str, data) -> str:
    """
    Save the extracted variant data in JSON format for the purposes of this codeathon.
    """

    # Save the data to a JSON file
    json_filename = gene + "_variants_extracted.json"
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)
    rprint(f"Saved variant information to {json_filename}")

    return json_filename


async def save_unique_conditions(gene: str, data) -> str:
    """
    Save the unique set of conditions for the given gene as a JSON file
    to be used downstream.
    """

    # Extract unique disease names
    unique_disease_names = {item["condition"]["text"] for item in data}

    # Save unique disease names to JSON
    unique_disease_json_filename = gene + "_unique_diseases.json"

    with open(unique_disease_json_filename, "w", encoding="utf-8") as json_file:
        # Convert the set to a list for JSON serialization
        json.dump(list(unique_disease_names), json_file, ensure_ascii=False, indent=4)
    print(f"Saved unique disease names to {unique_disease_json_filename}")

    # Assigning an ID to each unique condition and formatting the data
    formatted_conditions = [
        {"cid": idx + 1, "condition": condition_text}
        for idx, condition_text in enumerate(unique_disease_names)
    ]

    # Saving the formatted unique condition texts to a JSON file
    json_filename = gene + "_formatted_unique_conditions.json"
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(formatted_conditions, json_file, ensure_ascii=False, indent=4)

    print(f"Saved formatted unique condition texts to {json_filename}")

    return json_filename


async def main() -> None:
    """
    Coordinate the flow of data through the above functions within an
    asynchronous runtime that separates network/IO-bound tasks from
    CPU-bound tasks.
    """

    # Example usage parameters
    gene_file, gene, email = parse_command_line_args()

    if gene is not None:
        extracted_data = await pull_and_extract_data(gene, email)
        # _ = await save_variant_xml(gene, extracted_data)
        _ = await save_variant_json(gene, extracted_data)
        _ = await save_unique_conditions(gene, extracted_data)
        return

    # collect the list of genes
    with open(gene_file, "r", encoding="utf8") as input_handle:
        genes = [gene.strip() for gene in input_handle]

    # process each gene asynchronously
    for gene in genes:
        extracted_data = await pull_and_extract_data(gene, email)
        # _ = await save_variant_xml(gene, extracted_data)
        _ = await save_variant_json(gene, extracted_data)
        _ = await save_unique_conditions(gene, extracted_data)
        rprint(f"Data retrieval, extraction, and writing complete for {gene}")


if __name__ == "__main__":
    asyncio.run(main())
