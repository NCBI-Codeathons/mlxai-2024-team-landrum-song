import json


def load_json_data(file_path):
    """Load JSON data from a file."""
    with open(file_path, "r") as file:
        return json.load(file)


def update_variants_with_clusters(variants, clusters, print_excluded_RCVs):
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
    # print(len(exclude_from_clustering))
    # print("Excluded RCVs due to multiple conditions")

    # Name each cluster by the first name
    cluster_name = {}

    for cluster in clusters:
        if cluster["items"]:  # Check if there are any items in the cluster
            first_condition_name = cluster["items"][0]["content"]
            cluster_name[cluster["id"]] = first_condition_name

    # Print the cluster_name dictionary to verify
    # for cluster_id, name in cluster_name.items():
    # print(f"Cluster ID: {cluster_id} has the first condition name: {name}")

    # Mapping from condition text to cluster ID where cluster ID >= 0
    condition_to_cluster = {}
    for cluster in clusters:
        if int(cluster["id"]) >= 0:
            for item in cluster["items"]:
                condition_to_cluster[item["content"]] = cluster["id"]
    # print(list(condition_to_cluster.items())[:5])

    rcvs_before = len(set([v["RCV"] for v in variants]))

    # Step 1: Gather unique variant IDs
    unique_variant_ids = set(v["variant_id"] for v in variants_data)

    # Step 2: For each unique variant ID, process the RCVs and conditions
    updated_variants = []
    for variant_id in unique_variant_ids:
        # Filter variants by the current variant_id
        variants_with_id = [v for v in variants_data if v["variant_id"] == variant_id]

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
    rcvs_after = len(set([v["RCV"] for v in variants]))

    return variants, rcvs_before, rcvs_after


# Example gene names - replace with your actual list of genes
gene_names = ["LDLR", "SCN5A", "TSC1", "KCNQ1", "USH2A"]

for gene_name in gene_names:
    # Load the JSON data for variants and clusters
    variants_file_path = f"{gene_name}_variants_extracted.json"
    clusters_file_path = f"{gene_name}_clusters.json"
    curated_clusters_file_path = f"{gene_name}_clusters_curated.json"

    variants_data = load_json_data(variants_file_path)
    clusters_data = load_json_data(clusters_file_path)
    clusters_data_curated = load_json_data(curated_clusters_file_path)

    # Update the variants based on the clusters, respecting the new requirements
    print_excluded_RCVs = True
    updated_variants, rcvs_before, rcvs_after = update_variants_with_clusters(
        variants_data, clusters_data, print_excluded_RCVs
    )
    # print([cluster['id'] for cluster in clusters_data_curated])

    variants_data = load_json_data(variants_file_path)
    clusters_data = load_json_data(clusters_file_path)
    clusters_data_curated = load_json_data(curated_clusters_file_path)
    print_excluded_RCVs = False
    (
        updated_variants_curated,
        rcvs_before_curated,
        rcvs_after_curated,
    ) = update_variants_with_clusters(
        variants_data, clusters_data_curated, print_excluded_RCVs
    )

    # Save the updated data with a new file name
    with open(f"{gene_name}_variants_extracted_aggregated_RCVs.json", "w") as json_file:
        json.dump(updated_variants, json_file, indent=4)

    # Save the RCV count summary
    with open(f"{gene_name}_RCV_counts.txt", "w") as count_file:
        count_file.write(f"RCVs before: {rcvs_before}\nRCVs after: {rcvs_after}")

    # Print out the RCV counts for confirmation
    print(f"{gene_name}: RCVs before: {rcvs_before}, RCVs after: {rcvs_after}")

    # Repeat for curated
    # Save the updated data with a new file name
    with open(
        f"{gene_name}_variants_extracted_aggregated_RCVs_curated.json", "w"
    ) as json_file:
        json.dump(updated_variants_curated, json_file, indent=4)

    # Save the RCV count summary
    with open(f"{gene_name}_RCV_counts.txt", "w") as count_file:
        count_file.write(
            f"RCVs before: {rcvs_before_curated}\nRCVs after: {rcvs_after_curated}"
        )

    # Print out the RCV counts for confirmation
    print(
        f"Curated {gene_name}: RCVs before: {rcvs_before_curated}, RCVs after: {rcvs_after_curated}"
    )
    print("\n")
