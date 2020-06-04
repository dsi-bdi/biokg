<img src="https://i.ibb.co/MkVQGCG/biokg-logo-blue.png" alt="biokg-logo" width="300"/>

A knowledge graph for relational learning on biological data.

<br>

# Compiling BioKG with docker

## Building docker image
```bash
sudo docker build . -t dsi-bdi/biokg
```
## Running the image
```bash
sudo docker run --rm -v <data_path>:/biokg/data -e DB_USER='<drugbank_username>' -e DB_PASS='<drugbank_password>' dsi-bdi/biokg:latest
```
- where <data_path> is the fully qualified path to your data folder

## Data Sources
The biokg is built using the following data sources.

| Source Database                                    | License Type | URL                                                                  |
|----------------------------------------------------|--------------|----------------------------------------------------------------------|
| [UniProt](https://www.uniprot.org)                 | CC BY 4.0    | https://www.uniprot.org/help/license                                 |
| [Drugbank](https://www.drugbank.ca/)               | CC BY NC 4.0 | https://www.drugbank.ca/legal/terms_of_use                           |
| [KEGG](https://www.genome.jp/kegg/)                | Custom       | https://www.kegg.jp/kegg/legal.html                                  |
| [Sider](http://sideeffects.embl.de/)               | CC BY-NC-SA  | http://sideeffects.embl.de/about/                                    |
| [HPA](https://www.proteinatlas.org/)               | CC BY SA 3.0 | https://www.proteinatlas.org/about/licence                           |
| [Cellosaurus](https://web.expasy.org/cellosaurus/) | CC BY 4.0    | https://web.expasy.org/cgi-bin/cellosaurus/faq#Q22                   |
| [Reactome](https://reactome.org/)                  | CC0          | https://reactome.org/license                                         |
| [CTD](http://ctdbase.org/)                         | Custom       | http://ctdbase.org/about/legal.jsp                                   |
| [Intact](https://www.ebi.ac.uk/intact/)            | Apache 2.0   | https://www.ebi.ac.uk/intact/downloads                               |
| [MedGen](https://www.ncbi.nlm.nih.gov/medgen)      | Custom       | https://www.nlm.nih.gov/databases/download/terms_and_conditions.html |
| [MESH](https://www.ncbi.nlm.nih.gov/mesh)          | Custom       | https://www.nlm.nih.gov/databases/download/terms_and_conditions.html |
| [InterPro](http://www.ebi.ac.uk/interpro/)         | Custom       | ftp://ftp.ebi.ac.uk/pub/databases/interpro/release_notes.txt         |

# Funding
The development of this module has been fully supported by the CLARIFY project that has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 875160.
