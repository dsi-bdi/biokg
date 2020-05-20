from os import listdir
from os.path import isfile, join
from collections import defaultdict
import pandas as pd
from tabulate import tabulate


def summarize_folder(data_dp, data_name):
    print('##################################################################')
    print(f'#                         {data_name}                            ')
    print('##################################################################')
    data_files = [f for f in listdir(data_dp) if isfile(join(data_dp, f))]
    print(tabulate({"Files": data_files}, headers="keys"))
    print()


def summarize_file(filename, fp, columns, groups=[]):
    df = pd.read_csv(fp, sep='\t', names=columns)

    print(filename)
    print('\t'.join(f'<{c}>' for c in columns))

    for groupname in groups:
        group = df.groupby(groupname).size().to_dict()
        grouptitle = groupname
        if type(groupname) is list:
            grouptitle = ' - '.join(groupname)
        types = list(group.keys())
        types.append('TOTAL')

        counts = list(group.values())
        counts.append(len(df))
        print(tabulate(
            {
                grouptitle: types,
                'Counts': counts
            },
            headers="keys", tablefmt='orgtbl')
        )
        print()

    if len(groups) == 0:
        types = ['TOTAL']
        counts = [len(df)]
        print(tabulate(
            {
                'Triples': types,
                'Counts': counts
            },
            headers="keys", tablefmt='orgtbl')
        )
        print()


def uniprot_summary():
    summarize_folder(join(data_root, 'uniprot'), 'UNIPROT')
    print('Format:')
    summarize_file(
        'uniprot_metadata.txt',
        join(data_root, 'uniprot', 'uniprot_metadata.txt'),
        ['accession', 'metadata_type', 'value'],
        ['metadata_type']
    )

    summarize_file(
        'uniprot_facts.txt',
        join(data_root, 'uniprot', 'uniprot_facts.txt'),
        ['accession', 'fact_type', 'value'],
        ['fact_type']
    )
    summarize_file(
        'uniprot_ppi.txt',
        join(data_root, 'uniprot', 'uniprot_ppi.txt'),
        ['src_accession', 'INTERACTS_WITH', 'dest_accession']
    )


def cellosaurus_summary():
    summarize_folder(join(data_root, 'cellosaurus'), 'CELLOSAURUS')
    print('Format:')

    summarize_file(
        'cl_cat.txt',
        join(data_root, 'cellosaurus', 'cl_cat.txt'),
        ['cellosaurus_accession', 'names', 'species', 'sex', 'type']
    )

    summarize_file(
        'cl_map.txt',
        join(data_root, 'cellosaurus', 'cl_map.txt'),
        ['name', 'cellosaurus_accession']
    )

    summarize_file(
        'cl_pmid.txt',
        join(data_root, 'cellosaurus', 'cl_pmid.txt'),
        ['cellosaurus_accession', 'pubmed_ids']
    )

    summarize_file(
        'cl_geo.txt',
        join(data_root, 'cellosaurus', 'cl_geo.txt'),
        ['cellosaurus_accession', 'geo_id']
    )


def ctd_summary():
    summarize_folder(join(data_root, 'ctd'), 'CTD')
    print('Format:')
    # DISEASE
    summarize_file(
        'ctd_disease_biological_process.txt',
        join(data_root, 'ctd', 'ctd_disease_biological_process.txt'),
        ['disease_id', 'BIOLOGICAL_PROCESS', 'go_term_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_disease_biological_process.txt',
        join(data_root, 'ctd', 'ctd_disease_biological_process.txt'),
        ['disease_id', 'CELLULAR_COMPONENT', 'go_term_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_disease_molecular_function.txt',
        join(data_root, 'ctd', 'ctd_disease_molecular_function.txt'),
        ['disease_id', 'MOLECULAR_FUNCTION', 'go_term_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_disease_names.txt',
        join(data_root, 'ctd', 'ctd_disease_names.txt'),
        ['disease_id', 'DISEASE_NAME', 'name']
    )
    summarize_file(
        'ctd_disease_kegg_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_disease_kegg_pathway_association.txt'),
        ['disease_id', 'ASSOCIATED_PATHWAY', 'kegg_pathway_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_disease_reactome_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_disease_reactome_pathway_association.txt'),
        ['disease_id', 'ASSOCIATED_PATHWAY', 'reactome_pathway_id', 'data_status'],
        ['data_status']
    )

    # DRUG
    summarize_file(
        'ctd_drug_kegg_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_drug_kegg_pathway_association.txt'),
        ['drug_id', 'ASSOCIATED_PATHWAY', 'kegg_pathway_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_drug_reactome_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_drug_reactome_pathway_association.txt'),
        ['drug_id', 'ASSOCIATED_PATHWAY', 'reactome_pathway_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_drug_phenotype.txt',
        join(data_root, 'ctd', 'ctd_drug_phenotype.txt'),
        ['drug_id', 'DRUG_PHENOTYPE', 'go_term_id', 'effect', 'data_status', 'pubmed_ids'],
        ['effect', 'data_status']
    )
    summarize_file(
        'ctd_drug_protein_interactions.txt',
        join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt'),
        ['drug_id', 'interaction_type', 'accession', 'data_status', 'pubmed_ids'],
        ['interaction_type', 'data_status']
    )
    summarize_file(
        'ctd_drug_diseases.txt',
        join(data_root, 'ctd', 'ctd_drug_diseases.txt'),
        ['drug_id', 'ASSOCIATED_DISEASE', 'disease', 'data_status', 'pubmed_ids'],
        ['data_status']
    )

    # PROTEIN
    summarize_file(
        'ctd_protein_disease_association.txt',
        join(data_root, 'ctd', 'ctd_protein_disease_association.txt'),
        ['accession', 'ASSOCIATED_DISEASE', 'disease_id', 'data_status', 'pubmed_ids'],
        ['data_status']
    )
    summarize_file(
        'ctd_protein_kegg_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_protein_kegg_pathway_association.txt'),
        ['accession', 'ASSOCIATED_PATHWAY', 'reactome_pathway_id', 'data_status'],
        ['data_status']
    )
    summarize_file(
        'ctd_protein_reactome_pathway_association.txt',
        join(data_root, 'ctd', 'ctd_protein_reactome_pathway_association.txt'),
        ['accession', 'ASSOCIATED_PATHWAY', 'reactome_pathway_id', 'data_status'],
        ['data_status']
    )


def drugbank_summary():
    summarize_folder(join(data_root, 'drugbank'), 'DRUGBANK')
    print('Format:')
    summarize_file(
        'db_meta.txt',
        join(data_root, 'drugbank', 'db_meta.txt'),
        ['drug_id', 'metadata_type', 'value'],
        ['metadata_type']
    )
    summarize_file(
        'db_atc.txt',
        join(data_root, 'drugbank', 'db_atc.txt'),
        ['drug_id', 'DRUG_ATC_C[1-5]', 'atc_code']
    )
    summarize_file(
        'db_classification.txt',
        join(data_root, 'drugbank', 'db_classification.txt'),
        ['drug_id', 'classification_type', 'classification'],
        ['classification_type']
    )
    summarize_file(
        'db_product_stage.txt',
        join(data_root, 'drugbank', 'db_product_stage.txt'),
        ['drug_id', 'PRODUCT_STAGE', 'stage'],
        ['stage']
    )
    summarize_file(
        'db_mesh.txt',
        join(data_root, 'drugbank', 'db_mesh.txt'),
        ['drug_id', 'MESH_CATEGORY', 'mesh_id']
    )
    summarize_file(
        'db_mechanism_of_action.txt',
        join(data_root, 'drugbank', 'db_mechanism_of_action.txt'),
        ['drug_id', 'mechanism of action string']
    )
    summarize_file(
        'db_pathways.txt',
        join(data_root, 'drugbank', 'db_pathways.txt'),
        ['drug_id', 'predicate', 'object'],
        ['predicate']
    )
    summarize_file(
        'db_targets.txt',
        join(data_root, 'drugbank', 'db_targets.txt'),
        ['drug_id', 'drug_target', 'accession', 'action', 'pubmed_ids'],
        ['drug_target', 'action']
    )
    summarize_file(
        'db_interactions.txt',
        join(data_root, 'drugbank', 'db_interactions.txt'),
        ['src_drug_id', 'DRUG_INTERACTION', 'dest_drug_id', 'effect'],
        ['effect']
    )


def hpa_summary():
    summarize_folder(join(data_root, 'hpa'), 'HPA')
    print('Format:')

    summarize_file(
        'hpa_antibodies.txt',
        join(data_root, 'hpa', 'hpa_antibodies.txt'),
        ['antibody_id', 'accession', 'antigen_sequence']
    )

    summarize_file(
        'hpa_tissues_exp.txt',
        join(data_root, 'hpa', 'hpa_tissues_exp.txt'),
        ['accssion', 'tissue', 'expression level']
    )

    summarize_file(
        'hpa_cellines_exp.txt',
        join(data_root, 'hpa', 'hpa_cellines_exp.txt'),
        ['accssion', 'celline', 'organ#cellosaurus_accesion', 'expression level']
    )


def intact_summary():
    summarize_folder(join(data_root, 'intact'), 'INTACT')
    print('Format:')

    summarize_file(
        'intact_ppi.txt',
        join(data_root, 'intact', 'intact_ppi.txt'),
        ['src_accession', 'INTERACTS_WITH', 'dest_accession', 'pubmed_ids']
    )


def kegg_summary():
    summarize_folder(join(data_root, 'kegg'), 'KEGG')
    print('Format:')
    summarize_file(
        'disease_pathway.txt',
        join(data_root, 'kegg', 'disease_pathway.txt'),
        ['disease', 'DISEASE_PATHWAY', 'pathway']
    )

    summarize_file(
        'disease_drug.txt',
        join(data_root, 'kegg', 'disease_drug.txt'),
        ['disease', 'DISEASE_DRUG', 'drug']
    )

    summarize_file(
        'glycan_pathway.txt',
        join(data_root, 'kegg', 'glycan_pathway.txt'),
        ['glycan', 'GLYCAN_PATHWAY', 'pathway']
    )

    summarize_file(
        'drug_pathway.txt',
        join(data_root, 'kegg', 'drug_pathway.txt'),
        ['drug', 'DRUG_PATHWAY', 'pathway']
    )

    summarize_file(
        'network_pathway.txt',
        join(data_root, 'kegg', 'network_pathway.txt'),
        ['network', 'NETWORK_PATHWAY', 'pathway']
    )

    summarize_file(
        'network_drug.txt',
        join(data_root, 'kegg', 'network_drug.txt'),
        ['network', 'NETWORK_DRUG', 'drug']
    )

    summarize_file(
        'network_disease.txt',
        join(data_root, 'kegg', 'network_disease.txt'),
        ['network', 'NETWORK_DDISEASE', 'disease']
    )

    summarize_file(
        'gene_pathway.txt',
        join(data_root, 'kegg', 'gene_pathway.txt'),
        ['gene', 'GENE_PATHWAY', 'pathway']
    )

    summarize_file(
        'gene_disease.txt',
        join(data_root, 'kegg', 'gene_disease.txt'),
        ['gene', 'GENE_DISEASE', 'disease']
    )

    summarize_file(
        'gene_drug.txt',
        join(data_root, 'kegg', 'gene_drug.txt'),
        ['gene', 'GENE_DRUG', 'drug']
    )

    summarize_file(
        'gene_network.txt',
        join(data_root, 'kegg', 'gene_network.txt'),
        ['gene', 'GENE_NETWORK', 'network']
    )


def phosphosite_summary():
    summarize_folder(join(data_root, 'phosphosite'), 'PHOSPHOSITE')
    print('Format:')
    summarize_file(
        'phosphorylation_site.txt',
        join(data_root, 'phosphosite', 'phosphorylation_site.txt'),
        ['accession', 'PHOSPHORYLATION_SITE', 'site', 'species'],
        ['species']
    )

    summarize_file(
        'kinase_substrate.txt',
        join(data_root, 'phosphosite', 'kinase_substrate.txt'),
        ['kinase accession', 'PHOSPHORYLATES', 'substrate accession', 'site', 'kinase species', 'substrate species'],
        [['kinase species', 'substrate species']]
    )


def reactome_summary():
    summarize_folder(join(data_root, 'reactome'), 'REACTOME')


def sider_summary():
    summarize_folder(join(data_root, 'sider'), 'SIDER')
    print('Format:')
    summarize_file(
        'sider_effects_meta.txt',
        join(data_root, 'sider', 'sider_effects_meta.txt'),
        ['sider_effect_id', 'NAME', 'name string']
    )
    summarize_file(
        'sider_effects.txt',
        join(data_root, 'sider', 'sider_effects.txt'),
        ['sider_id', 'SIDE_EFFECT', 'sider_effect_id']
    )
    summarize_file(
        'sider_indications_meta.txt',
        join(data_root, 'sider', 'sider_indications_meta.txt'),
        ['sider_indication_id', 'NAME', 'name string']
    )
    summarize_file(
        'sider_indications.txt',
        join(data_root, 'sider', 'sider_indications.txt'),
        ['sider_id', 'SIDE_EFFECT', 'sider_indication_id']
    )


data_root = 'data/preprocessed'


if __name__ == '__main__':
    uniprot_summary()
    ctd_summary()
    drugbank_summary()
    kegg_summary()
    phosphosite_summary()
    reactome_summary()
    sider_summary()
    intact_summary()
    hpa_summary()
    cellosaurus_summary()
