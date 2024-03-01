#!/usr/bin/env python3

"""
```
usage: aggregate_by_rcv.py [-h] --variants_path VARIANTS_PATH --cluster_path CLUSTER_PATH --gene GENE

options:
  -h, --help            show this help message and exit
  --variants_path VARIANTS_PATH, -v VARIANTS_PATH
                        Path to JSON of all RCV variants alongside the associated conditions.
  --cluster_path CLUSTER_PATH, -c CLUSTER_PATH
                        Path to a JSON of condition clusters produced with llm
  --gene GENE, -g GENE  Gene of interest
```
"""

import argparse
import json

from pathlib import Path
from typing import Optional, Tuple


def parse_command_line_args() -> Tuple[Path, Path, Optional[Path], str]:
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--variants_path",
        "-v",
        type=Path,
        required=True,
        help="Path to JSON of all RCV variants alongside the associated conditions.",
    )
    parser.add_argument(
        "--cluster_path",
        "-c",
        type=Path,
        required=True,
        help="Path to a JSON of condition clusters produced with llm",
    )
    parser.add_argument(
        "--curated_path",
        "-C",
        type=Optional[Path],
        required=False,
        default=None,
        help="Path to a JSON of condition clusters produced with llm",
    )
    parser.add_argument(
        "--gene",
        "-g",
        type=str,
        required=True,
        help="Gene of interest",
    )
    args = parser.parse_args()

    return args.variants_path, args.cluster_path, args.curated_path, args.gene


def load_json_data(file_path):
    """Load JSON data from a file."""
    with open(file_path, "r", encoding="utf8") as file:
        return json.load(file)


def update_variants_with_clusters(
    variants, clusters, print_excluded_RCVs: bool, gene_name: str
):
    """ """

    # Identify RCVs with multiple different conditions (for exclusion from clustering, but not from counts)
    rcv_condition_counts = {}
    for variant in variants:
        rcv = variant["RCV"]
        if rcv in rcv_condition_counts:
            rcv_condition_counts[rcv].add(variant["condition"]["text"])
        else:
            rcv_condition_counts[rcv] = {variant["condition"]["text"]}

    # RCVs to exclude from clustering due to multiple conditions
    exclude_from_clustering = {
        rcv for rcv, conditions in rcv_condition_counts.items() if len(conditions) > 1
    }
    if print_excluded_RCVs:
        print(
            f"Number of excluded RCVs: {len(exclude_from_clustering)} for {gene_name}"
        )
        total_conditions_in_excluded_rcvs = sum(
            len(conditions)
            for rcv, conditions in rcv_condition_counts.items()
            if rcv in exclude_from_clustering
        )
        print(
            f"Total number of conditions in excluded RCVs: {total_conditions_in_excluded_rcvs} for {gene_name}"
        )

    # Name each cluster by the first name
    cluster_name = {}

    for cluster in clusters:
        if cluster["items"]:  # Check if there are any items in the cluster
            first_condition_name = cluster["items"][0]["content"]
            cluster_name[cluster["id"]] = first_condition_name

    # Mapping from condition text to cluster ID where cluster ID >= 0
    condition_to_cluster = {}
    for cluster in clusters:
        if int(cluster["id"]) >= 0:
            for item in cluster["items"]:
                condition_to_cluster[item["content"]] = cluster["id"]

    rcvs_before = len(set(v["RCV"] for v in variants))

    # Step 1: Gather unique variant IDs
    unique_variant_ids = set(v["variant_id"] for v in variants)

    # Step 2: For each unique variant ID, process the RCVs and conditions
    updated_variants = []
    for variant_id in unique_variant_ids:
        # Filter variants by the current variant_id
        variants_with_id = [v for v in variants if v["variant_id"] == variant_id]

        # Mapping from condition to RCV for this variant_id
        condition_rcv_mapping = {
            v["condition"]["text"]: v["RCV"] for v in variants_with_id
        }

        # Check if multiple RCVs are part of the same cluster and update
        cluster_to_rcvs = {}
        for condition, rcv in condition_rcv_mapping.items():
            cluster_id = condition_to_cluster.get(condition, None)
            if cluster_id is not None:
                if cluster_id in cluster_to_rcvs:
                    cluster_to_rcvs[cluster_id].append(rcv)
                else:
                    cluster_to_rcvs[cluster_id] = [rcv]

        # If multiple RCVs are in the same cluster, choose one RCV to represent all
        for rcvs in cluster_to_rcvs.values():
            if len(rcvs) > 1:
                chosen_rcv = rcvs[0]  # Choose the first RCV
                for condition, rcv in condition_rcv_mapping.items():
                    if rcv in rcvs:
                        condition_rcv_mapping[condition] = chosen_rcv

        # Update the variants list
        for variant in variants_with_id:
            variant["RCV"] = condition_rcv_mapping[variant["condition"]["text"]]
            updated_variants.append(variant)

    # At this point, updated_variants contains the updated list of variants
    # Proceed to save this updated list to a JSON file or use as needed

    # Count unique RCVs before and after (including those with multiple conditions)
    rcvs_after = len(set(v["RCV"] for v in variants))

    return variants, rcvs_before, rcvs_after


def compare_with_curated(variants_data, clusters_data_curated, gene_name):
    print_excluded_RCVs = False
    (
        updated_variants_curated,
        rcvs_before_curated,
        rcvs_after_curated,
    ) = update_variants_with_clusters(
        variants_data, clusters_data_curated, print_excluded_RCVs, gene_name
    )

    # Save the updated data with a new file name
    with open(
        f"{gene_name}_variants_extracted_aggregated_RCVs_curated.json",
        "w",
        encoding="utf8",
    ) as json_file:
        json.dump(updated_variants_curated, json_file, indent=4)

    # Save the RCV count summary
    with open(f"{gene_name}_RCV_counts.txt", "w", encoding="utf8") as count_file:
        count_file.write(
            f"RCVs before: {rcvs_before_curated}\nRCVs after: {rcvs_after_curated}"
        )

    # Print out the RCV counts for confirmation
    print(
        f"Curated {gene_name}: RCVs before: {rcvs_before_curated}, RCVs after: {rcvs_after_curated}"
    )
    print("\n")


def main() -> None:
    """ """

    (
        variants_file_path,
        clusters_file_path,
        clusters_data_curated,
        gene_name,
    ) = parse_command_line_args()
    curated_clusters_file_path = clusters_file_path

    variants_data = load_json_data(variants_file_path)
    clusters_data = load_json_data(clusters_file_path)

    # Update the variants based on the clusters, respecting the new requirements
    print_excluded_RCVs = True
    updated_variants, rcvs_before, rcvs_after = update_variants_with_clusters(
        variants_data, clusters_data, print_excluded_RCVs, gene_name
    )

    if clusters_data_curated is not None:
        clusters_data_curated = load_json_data(curated_clusters_file_path)
        compare_with_curated(variants_data, clusters_data_curated, gene_name)

    # Save the updated data with a new file name
    with open(
        f"{gene_name}_variants_extracted_aggregated_RCVs.json", "w", encoding="utf8"
    ) as json_file:
        json.dump(updated_variants, json_file, indent=4)

    # Save the RCV count summary
    with open(f"{gene_name}_RCV_counts.txt", "w", encoding="utf8") as count_file:
        count_file.write(f"RCVs before: {rcvs_before}\nRCVs after: {rcvs_after}")

    # Print out the RCV counts for confirmation
    print(f"{gene_name}: RCVs before: {rcvs_before}, RCVs after: {rcvs_after}")


if __name__ == "__main__":
    main()
