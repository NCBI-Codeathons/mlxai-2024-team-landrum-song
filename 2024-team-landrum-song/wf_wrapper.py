#!/usr/bin/env python3

"""
TODO
"""

from flytekit import workflow, task
from rich import print as rprint


@task
def retrieve_data(query_gene: str) -> None:
    """
    Use the entrez Python API to retrieve data from NCBI ClinVar
    for a given gene. This is the chief source of network I/O in
    the workflow.
    """
    rprint(f"retrieving ClinVar data for {query_gene}")


@task
def extract_condition_set():
    """
    Parse the raw XML output from Entrez into a set of unique condition
    names whose structure may be variously nested.
    """
    rprint("Extracting set of conditions from the ClinVar query.")


@task
def prompt_llm():
    """
    Automatically synthesize a prompt for an LLM and send it in.
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


@workflow
def main() -> None:
    """
    TODO
    """


if __name__ == "__main__.py":
    main()
