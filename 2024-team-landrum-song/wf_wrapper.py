#!/usr/bin/env python3

"""
TODO
"""

import argparse
from pathlib import Path
from typing import Tuple

from flytekit import task, workflow
from rich import print as rprint


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
    args = parser.parse_args()

    return args.genelist, args.email


@task
def retrieve_data(query_gene: str, _email: str) -> None:
    """
    Use the entrez Python API to retrieve data from NCBI ClinVar
    for a given gene. This is the chief source of network I/O in
    the workflow.
    """
    rprint(f"retrieving ClinVar data for {query_gene}")

    # call wrapper function from entrez-xml script
    # parsed_xml = search_and_save_gene_info_combined(query_gene, email)


@task
def extract_condition_set() -> None:
    """
    Parse the raw XML output from Entrez into a set of unique condition
    names whose structure may be variously nested.
    """
    rprint("Extracting set of conditions from the ClinVar query.")


@task
def prompt_llm() -> None:
    """
    Automatically synthesize a prompt for an LLM and send it in,
    collecting its responses
    """
    condition_set = set(["test1", "test2"])

    prompt_template = """
    Can you take the following disease names and sort them into a nested
    JSON where diseases that are subtypes of other diseases in the list
    are nested within those diseases? I want to create a JSON map of
    disease types and subtypes so that diseases with slightly different
    names aren't treated as entirely different diseases in my database.
    """

    _llm_prompt = f"{prompt_template}\n\n{condition_set}"


@task
def parse_response() -> None:
    """
    TODO
    """


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
    genelist, email = parse_command_line_args()

    # run the raw data retrieval tasks
    retrieved_datasets = [retrieve_data(gene, email) for gene in genelist]

    # extract unique sets of conditions for each gene dataset
    condition_sets = [extract_condition_set() for dataset in retrieved_datasets]

    # send in the prompts
    llm_responses = [prompt_llm() for diseases in condition_sets]

    # parse responses into JSON formats
    sorted_conditions = [parse_response() for response in llm_responses]

    # remap the new top-level disease names onto the original clinvar data
    final_mappings = [map_back_to_clinvar() for cond in sorted_conditions]

    # write out results for manual review
    _ = [write_out_results() for mapping in final_mappings]


if __name__ == "__main__.py":
    main()
