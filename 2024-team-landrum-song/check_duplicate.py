#!/usr/bin/env python3

import sys
import json


GENE = sys.argv[1]


def deduplicate_in_clusters(gene: str):
    """
    TODO
    """

    # Load LDLR_test.json
    with open(gene + "_variants_extracted.json", "r", encoding="utf8") as file:
        ldlr_test_data = json.load(file)

    # Load LDLR_clusters.json
    with open(gene + "_clusters.json", "r", encoding="utf8") as file:
        ldlr_clusters_data = json.load(file)

    # Step 2: Organize LDLR_clusters_data for efficient searching
    # Exclude clusters with id: "-1" and create a mapping of content to cluster id

    content_to_cluster_id = {}
    for cluster in ldlr_clusters_data:
        if cluster["id"] != "-1":
            for item in cluster["items"]:
                content_to_cluster_id[item["content"].lower()] = cluster["id"]

    # Step 3: Process LDLR_test_data to find records meeting the criteria
    # For each variant_id, we will check if 2 or more text fields can be found in any content fields of clusters

    # Organize LDLR_test_data by variant_id and collect texts
    variant_id_to_texts = {}
    for record in ldlr_test_data:
        variant_id = record["variant_id"]
        text = record["condition"]["text"].lower()
        if variant_id not in variant_id_to_texts:
            variant_id_to_texts[variant_id] = []
        variant_id_to_texts[variant_id].append(text)

    # Find qualifying records
    qualifying_records = []
    total_reduced_records = 0
    counter_unique_variant_id_duplicate = (
        0  # count the number of unique variant_id with similar disease names
    )
    for variant_id, texts in variant_id_to_texts.items():
        # Count texts found in any content fields of clusters
        found_texts_count = sum(1 for text in texts if text in content_to_cluster_id)

        # If 2 or more texts found, add all records with this variant_id to qualifying_records
        if found_texts_count >= 2:
            qualifying_records.extend(
                [
                    record
                    for record in ldlr_test_data
                    if record["variant_id"] == variant_id
                ]
            )
            total_reduced_records += found_texts_count - 1
            counter_unique_variant_id_duplicate = (
                counter_unique_variant_id_duplicate + 1
            )

    print(
        "\nNumber of RCV records that can be potentially reduced: ",
        total_reduced_records,
    )
    print("Total nmber of RCV records: ", len(ldlr_test_data))

    print(
        "\nNumber of variants with groupable disease names: ",
        counter_unique_variant_id_duplicate,
    )
    print(
        "Number of unique variants: ", len(variant_id_to_texts)
    )  # number of unique variant_id
    print()

    # Step 4: Save the qualifying records into a new JSON file

    output_file_path = gene + "_qualifying_records.json"
    with open(output_file_path, "w", encoding="utf8") as outfile:
        json.dump(qualifying_records, outfile, indent=4)

    print("Write ", output_file_path)


if __name__ == "__main__":
    deduplicate_in_clusters(GENE)
