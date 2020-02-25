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
    swissprot_hsa_id_map_fp = join(sources_dp, "swissprot_human_ids_mapping.txt.gz")
    swissprot_hsa_ids_fp = join(sources_dp, "swissprot_human_ids.txt.gz")
    swissprot_hsa_entries_fp = join(sources_dp, "swissprot_human_entries.txt.gz")

    download_file_md5_check(srcs_cp["uniprot"]["swissprot_human_ids_mapping"], swissprot_hsa_id_map_fp)
    download_file_md5_check(srcs_cp["uniprot"]["swissprot_human_ids"], swissprot_hsa_ids_fp)
    download_file_md5_check(srcs_cp["uniprot"]["swissprot_human_entries"], swissprot_hsa_entries_fp)
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
    """ Download cellosaurus database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    print_section_header("Downloading Reactome data files")

    reactome_ppi_fp = join(sources_dp, "reactome_ppi.txt.gz")
    reactome_pathway_rels_fp = join(sources_dp, "reactome_pathway_rels.txt.gz")
    reactome_protein_complex_rels_fp = join(sources_dp, "reactome_protein_complex_rels.txt.gz")
    reactome_complex_pathway_rels_fp = join(sources_dp, "reactome_complex_pathway_rels.txt.gz")
    reactome_go_mapping_fp = join(sources_dp, "reactome_go_mapping.txt.gz")

    download_file_md5_check(srcs_cp["reactome"]["reactome_ppi"], reactome_ppi_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_pathway_rels"], reactome_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_protein_complex_rels"], reactome_protein_complex_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_complex_pathway_rels"], reactome_complex_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_go_mapping"], reactome_go_mapping_fp)

    print_bold_line()
