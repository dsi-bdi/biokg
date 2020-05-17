# -*- coding: utf-8 -*-

from configparser import RawConfigParser
from .util.io import *
from .util.extras import *
from os.path import join


def download_uniprot_files(sources_dp, srcs_cp):
    """ Download uniprot database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    # download source data (uniprot)
    print_section_header("Downloading UniProt data files")
    swissprot_entries_fp = join(sources_dp, "swissprot_entries.txt.gz")

    download_file_md5_check(srcs_cp["uniprot"]["swissprot_entries"], swissprot_entries_fp)
    print_bold_line()


def download_cellosaurus_data(sources_dp, srcs_cp):
    """ Download cellosaurus database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Cellosaurus data files")
    cellosaurus_data_fp = join(sources_dp, "cellosaurus_data.txt")
    download_file_md5_check(srcs_cp["cellosaurus"]["cellosaurus_txt"], cellosaurus_data_fp)
    print_bold_line()


def download_hpa_data(sources_dp, srcs_cp):
    """ Download cellosaurus database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Human Protein Atlas data files")
    hpa_data_fp = join(sources_dp, "proteinatlas.xml.gz")
    normal_tissues_fp = join(sources_dp, "normal_tissues.tsv.gz")
    rna_cellines_fp = join(sources_dp, "rna_cellines.tsv.gz")
    pathology_fp = join(sources_dp, "pathology.tsv.gz")

    download_file_md5_check(srcs_cp["hpa"]["xml_database"], hpa_data_fp)
    download_file_md5_check(srcs_cp["hpa"]["normal_tissues"], normal_tissues_fp)
    download_file_md5_check(srcs_cp["hpa"]["rna_cellines"], rna_cellines_fp)
    download_file_md5_check(srcs_cp["hpa"]["pathology"], pathology_fp)

    print_bold_line()


def download_reactome_data(sources_dp, srcs_cp):
    """ Download reactome database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Reactome data files")

    reactome_ppi_fp = join(sources_dp, "reactome_ppi.txt")
    reactome_pathway_rels_fp = join(sources_dp, "reactome_pathway_rels.txt")
    reactome_protein_complex_rels_fp = join(sources_dp, "reactome_protein_complex_rels.txt")
    reactome_complex_pathway_rels_fp = join(sources_dp, "reactome_complex_pathway_rels.txt")
    reactome_go_mapping_fp = join(sources_dp, "reactome_go_mapping.txt.gz")
    reactome_omim_mapping_fp = join(sources_dp, "reactome_omim_mapping.txt")

    download_file_md5_check(srcs_cp["reactome"]["reactome_ppi"], reactome_ppi_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_pathway_rels"], reactome_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_protein_complex_rels"], reactome_protein_complex_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_complex_pathway_rels"], reactome_complex_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_go_mapping"], reactome_go_mapping_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_omim_mapping"], reactome_omim_mapping_fp)

    print_bold_line()


def download_ctd_data(sources_dp, srcs_cp):
    """ Download ctd database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading CTD data files")

    ctd_chemical_id_mapping = join(sources_dp,'CTD_chemicals.tsv.gz')
    ctd_gene_id_mapping =join(sources_dp,'CTD_genes.tsv.gz')
    ctd_chemical_gene_interactions = join(sources_dp,'CTD_chem_gene_ixns.tsv.gz')
    ctd_chemical_disase_association =join(sources_dp,'CTD_chemicals_diseases.tsv.gz')
    ctd_chemical_pathway_association = join(sources_dp,'CTD_chem_pathways_enriched.tsv.gz')
    ctd_gene_disease_association = join(sources_dp,'CTD_genes_diseases.tsv.gz')
    ctd_gene_pathway_association = join(sources_dp,'CTD_genes_pathways.tsv.gz')
    ctd_disease_pathway_assiciation = join(sources_dp,'CTD_diseases_pathways.tsv.gz')
    ctd_disease_molecular_function = join(sources_dp,'CTD_disease_molecular_function.tsv.gz')
    ctd_disease_cellular_component = join(sources_dp,'CTD_disease_cellular_component.tsv.gz')
    ctd_disease_biological_process = join(sources_dp,'CTD_disease_biological_process.tsv.gz')
    ctd_chemical_phenotype = join(sources_dp,'CTD_chemical_phenotype.tsv.gz')

    download_file_md5_check(srcs_cp["ctd"]["chemical_id_mapping"], ctd_chemical_id_mapping)
    download_file_md5_check(srcs_cp["ctd"]["gene_id_mapping"], ctd_gene_id_mapping)
    download_file_md5_check(srcs_cp["ctd"]["chemical_gene_interactions"], ctd_chemical_gene_interactions)
    download_file_md5_check(srcs_cp["ctd"]["chemical_disase_association"], ctd_chemical_disase_association)
    download_file_md5_check(srcs_cp["ctd"]["chemical_pathway_association"], ctd_chemical_pathway_association)
    download_file_md5_check(srcs_cp["ctd"]["gene_disease_association"], ctd_gene_disease_association)
    download_file_md5_check(srcs_cp["ctd"]["gene_pathway_association"], ctd_gene_pathway_association)
    download_file_md5_check(srcs_cp["ctd"]["disease_pathway_assiciation"], ctd_disease_pathway_assiciation)
    download_file_md5_check(srcs_cp["ctd"]["disease_molecular_function"], ctd_disease_molecular_function)
    download_file_md5_check(srcs_cp["ctd"]["disease_cellular_component"], ctd_disease_cellular_component)
    download_file_md5_check(srcs_cp["ctd"]["disease_biological_process"], ctd_disease_biological_process)
    download_file_md5_check(srcs_cp["ctd"]["chemical_phenotype"], ctd_chemical_phenotype)

    print_bold_line()
    
def download_phosphosite_data(sources_dp, srcs_cp):
    """ Download phosphositeplus database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading PhosphoSitePlus data files")

    phosphorylation_site_fp = join(sources_dp, "phosphorylation_site.txt.gz")
    kinase_substrate_fp = join(sources_dp, "kinase_substrate.txt.gz")

    download_file_md5_check(srcs_cp["phosphositeplus"]["phosphorylation_site"], phosphorylation_site_fp)
    download_file_md5_check(srcs_cp["phosphositeplus"]["kinase_substrate"], kinase_substrate_fp)

    print_bold_line()


def download_intact_data(sources_dp, srcs_cp):
    """ Download intact database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Intact data files")

    intact_zip_fp = join(sources_dp, "intact.zip")

    download_file_md5_check(srcs_cp["intact"]["intact_zip"], intact_zip_fp)

    print_bold_line()


def download_sider_data(sources_dp, srcs_cp):
    """ Download sider database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading SIDER data files")

    sider_interactions_fp = join(sources_dp, "sider_interactions.tsv.gz")
    sider_side_effects_fp = join(sources_dp, "sider_side_effects.tsv.gz")

    download_file_md5_check(srcs_cp["sider"]["indications"], sider_interactions_fp)
    download_file_md5_check(srcs_cp["sider"]["side_effects"], sider_side_effects_fp)

    print_bold_line()


def download_drugbank_data(sources_dp, srcs_cp, username, password):
    """ Download drugbank database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Drugbank data files")

    full_database_fp = join(sources_dp, "drugbank_all_full_database.xml.zip")
    download_file_md5_check(
        srcs_cp["drugbank"]["drugbank_all_full_database"], 
        full_database_fp, 
        username = username, 
        password = password
    )

    print_bold_line()