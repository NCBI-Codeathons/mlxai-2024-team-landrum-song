#!/usr/bin/env python3

"""
TODO
"""

import sys
import argparse
import functools
import subprocess
from pathlib import Path
from typing import Dict, List, Set, Tuple

from flytekit import map_task, task, workflow
from flytekit.types.file import FlyteFile
from rich.pretty import pprint as rprint

from .pull_and_extract_clinvar import pull_and_extract_data, save_variant_json
from .check_duplicate import deduplicate_in_clusters


def parse_command_line_args() -> Tuple[Path, str]:
    """
    TODO
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
    parser.add_argument(
        "--tuning-param",
        "-t",
        type=float,
        required=False,
        default=1.0,
        help="Email for NCBI ClinVar tracking purposes.",
    )
    args = parser.parse_args()

    return args.genelist, args.email


@task
def retrieve_data(query_gene: str, _email: str) -> Tuple[str, str]:
    """
    Use the entrez Python API to retrieve data from NCBI ClinVar
    for a given gene. This is the chief source of network I/O in
    the workflow.
    """
    rprint(f"retrieving ClinVar data for {query_gene}")

    # call wrapper function from entrez-xml script, and
    # extract the relevant data as a list of dictionaries
    extracted_data = pull_and_extract_data(query_gene, _email)

    # write to a JSON
    json_out = save_variant_json(query_gene, extracted_data)

    return query_gene, json_out


@task
def extract_condition_set(data: List[Dict]) -> Set[str]:
    """
    Parse the raw XML output from Entrez into a set of unique condition
    names whose structure may be variously nested.
    """
    rprint("Extracting set of conditions from the ClinVar query.")

    condition_set = set(record.get("condition").get("text") for record in data)

    return condition_set


@task
def prompt_llm(_condition_set: Set[str]) -> FlyteFile:
    """
    Prompt a language model to perform Density-Based Spatial
    Clustering of Applications with Noise (DBSCAN) on the condition
    set for the current gene and extract cluster indices.
    """

    cluster_json_path = "condition_clusters.json"

    # command 1
    command1 = (
        "llm",
        "embed-multi",
        "diseases",
        "test2.json",
        "--database",
        "tmp.db",
        "--model",
        "sentence-transformers/all-MiniLM-L6-v2",
        "--store",
    )
    llm_embed_process = subprocess.Popen(command1, stdout=subprocess.PIPE)

    # command 2
    command2 = ("llm", "cluster", "diseases", "--database", "tmp.db", "1")
    llm_cluster_process = subprocess.Popen(command2, stdout=cluster_json_path)

    # Allow vcftools to receive a SIGPIPE if gzip exits
    if llm_embed_process.stdout:
        llm_embed_process.stdout.close()

    # Wait for processes to complete and get their exit codes and stderr
    _, embed_stderr = llm_cluster_process.communicate()
    _, cluster_stderr = llm_cluster_process.communicate()

    if llm_embed_process.returncode != 0 or llm_cluster_process.returncode != 0:
        rprint("LLM clustering failed. See errors below:")
        rprint(f"Embedding stderr:\{embed_stderr}")
        rprint(f"Clustering stderr:\{cluster_stderr}")
        sys.exit(1)

    return FlyteFile(path=cluster_json_path)


@task
def parse_and_dedup_response(_cluster_json: FlyteFile) -> None:
    """
    TODO
    """

    deduplicate_in_clusters("test")


@task
def map_back_to_clinvar() -> None:
    """
    TODO
    """


@task
def write_out_results() -> None:
    """
    TODO
    """


@workflow
def main() -> None:
    """
    TODO
    """

    # parse command line arguments that supply a file listing genes of interest
    # and an email for NCBI tracking purposes
    gene_file, email = parse_command_line_args()

    # collect the list of genes
    with open(gene_file, "r", encoding="utf8") as input_handle:
        genes = [gene.strip() for gene in input_handle]

    # run the raw data retrieval tasks
    retrieval_partial = functools.partial(retrieve_data, email=email)
    retrieval_result = map_task(retrieval_partial)(query_gene=genes)

    # extract unique sets of conditions for each gene dataset
    condition_sets = map_task(extract_condition_set)(retrieval_result)

    # send in the prompts
    llm_responses = map_task(prompt_llm)(condition_sets)

    # parse responses into JSON formats
    sorted_conditions = map_task(parse_and_dedup_response)(llm_responses)

    # remap the new top-level disease names onto the original clinvar data
    _final_mappings = map_task(map_back_to_clinvar)(sorted_conditions)

    # write out results for manual review
    # _ = map_task(write_out_results)(final_mappings)


if __name__ == "__main__":
    main()
