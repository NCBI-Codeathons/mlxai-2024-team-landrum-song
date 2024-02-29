# ClinCluster: A Package for Aggregating Disease Terms in ClinVar

[![Docker CI](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/actions/workflows/docker-image.yaml)

List of participants and affiliations:
- Melissa Landrum, NIH/NCBI (Co-Team Leader)
- Guangfeng Song, NIH/NCBI (Co-Team Leader)
- Lauren Edgar, NIH/NHGRI
- Benjamin Kesler, Vanderbilt University
- Nicholas Minor, University of Wisconsin
- Michael Muchow, Unaffiliated
- Rebecca Orris, NIH/NCBI
- Wengang Zhang, NIH/NCI

## Project Goals
Naming for human genetic diseases is complex. Diseases may be named for the phenotype; other information may be alluded to including the relevant gene, mode of inheritance, or the mechanism of disease. Diseases may be described at a high level with a generic name, or at a lower level with a more specific name; however, in the context of a variant in a specific gene, these differences may not be considered important. ClinVar data would be easier to ingest in bulk and to read in web displays if there were a meaningful way to aggregate diseases that effectively mean the same thing in the context of a gene.
## Approach
Develop a ML/AI approach to aggregating diseases for 5 genes with variants in ClinVar: LDLR, KCNQ1, USH2A, SCN5A, TSC1.

* A proposed solution can be provided

* Use the gene symbol to guide whether different terms are meaningfully different for that gene, e.g. Familial hypercholesterolemia vs Hypercholesterolemia, familial, 1

* Use other information such as mode of inheritance, clinical features, mechanism of disease to decide if terms should be aggregated or not

If time, demonstrate how RCV records for variant-condition pairs in ClinVar for one or more genes would be different using aggregated diseases.

<img src="https://files.slack.com/files-pri/T06HWPPPTT8-F06M8R2JEE5/ncbi-clinvar-llm-dag.png" alt="NCBI-clinvar-llm-dag.png"/>![image](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/assets/34135674/bf3d92be-a132-4af9-bec1-87b287966d5b)

<img src="blob:chrome-untrusted://media-app/03eb0034-22d4-4894-8461-42ecad647455" alt="Use LLM to aggregate these similar disease terms.jpg"/>![image](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/assets/34135674/2eb07c2a-a678-40a9-b584-70073df4e426)

<img src="blob:chrome-untrusted://media-app/8428382f-fad0-428a-a345-3ed8c32cd871" alt="DBSCAN_ cluster the disease names.jpg"/>![image](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/assets/34135674/fcbfacd5-d0bc-4d14-8368-1e02b8e53c63)

### Aggregating the disease terms
1. Extract disease names from all RCV records associated with a single gene
  1. Extract all unique MedGen ID and their disease names

2. Cluster the name into their umbrella disease category
  1. Program LLM to cluster the disease terms
  2. Fork https://github.com/simonw/llm-cluster/blob/main/llm_cluster.py

3. Identify variant records with similar disease names that belong to the same umbrella diseases

4. Assign those RCV records with the corrected (umbrella) name

## Results

<img src="blob:chrome-untrusted://media-app/56f48609-672e-48c8-8364-1eb1882696e4" alt="Performance Metrics .jpg"/>![image](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/assets/34135674/b9c5d6d7-64a1-4f3d-90b2-577d2a49a7cc)


## Future Work

## Common Acronyms
Abbreviation  | Acronym
------------- | -------------
ACMG  | Content Cell
DBSCAN  | Content Cell
GPT  | Generative Pre-trained Transformer


Abbreviation | Acronym
------------   -------
ACMG         | American College of Medical Genetics and Genomics
DBSCAN       | Density-Based Spatial Clustering of Applications with Noise

## NCBI Codeathon Disclaimer
This software was created as part of an NCBI codeathon, a hackathon-style event focused on rapid innovation. While we encourage you to explore and adapt this code, please be aware that NCBI does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)
