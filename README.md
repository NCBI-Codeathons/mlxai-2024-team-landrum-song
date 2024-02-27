# Aggregation of disease terms in ClinVar

List of participants and affiliations:
- Melissa Landrum, NIH/NCBI (Co-Team Leader)
- Lauren Edgar, NIH/NHGRI
- Benjamin Kesler, Vanderbilt University
- Nicholas Minor, University of Wisconsin
- Michael Muchow, Unaffiliated
- Rebecca Orris, NIH/NCBI
- Hasitha Premathilake, NIH/NIA
- Guangfeng Song, NIH/NCBI
- Wengang Zhang, NIH

## Project Goals
Naming for human genetic diseases is complex. Diseases may be named for the phenotype; other information may be alluded to including the relevant gene, mode of inheritance, or the mechanism of disease. Diseases may be described at a high level with a generic name, or at a lower level with a more specific name; however, in the context of a variant in a specific gene, these differences may not be considered important. ClinVar data would be easier to ingest in bulk and to read in web displays if there were a meaningful way to aggregate diseases that effectively mean the same thing in the context of a gene. 
## Approach
Develop a ML/AI approach to aggregating diseases for 5 genes with variants in ClinVar: LDLR, KCNQ1, USH2A, SCN5A, TSC1. 

* A proposed solution can be provided

* Use the gene symbol to guide whether different terms are meaningfully different for that gene, e.g. Familial hypercholesterolemia vs Hypercholesterolemia, familial, 1

* Use other information such as mode of inheritance, clinical features, mechanism of disease to decide if terms should be aggregated or not 

If time, demonstrate how RCV records for variant-condition pairs in ClinVar for one or more genes would be different using aggregated diseases.

<img src="https://files.slack.com/files-pri/T06HWPPPTT8-F06M8R2JEE5/ncbi-clinvar-llm-dag.png" alt="NCBI-clinvar-llm-dag.png"/>![image](https://github.com/NCBI-Codeathons/mlxai-2024-team-landrum-song/assets/34135674/bf3d92be-a132-4af9-bec1-87b287966d5b)

## Results

## Future Work

## NCBI Codeathon Disclaimer
This software was created as part of an NCBI codeathon, a hackathon-style event focused on rapid innovation. While we encourage you to explore and adapt this code, please be aware that NCBI does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)

