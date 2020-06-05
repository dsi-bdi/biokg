# -*- coding: utf-8 -*-

from os.path import join, isdir, exists
from os import mkdir, environ
from configparser import RawConfigParser
import sys
from biokg.loader import *
from biokg.util.extras import program_header
from biokg.processing.parsers import *
from biokg.util.io import export_file_md5, file_has_valid_md5


def preprocess_graph():
    """ Program entry point
    """
    print_bold_line()
    print(program_header)
    print_bold_line()

    # define work directories and files
    data_dp = "./data/"
    sources_dp = join(data_dp, "sources")
    preprocessed_dp = join(data_dp, "preprocessed")
    output_dp = join(data_dp, "output")
    sources_fp = "./sources.ini"

    # create directories if not existing
    mkdir(data_dp) if not isdir(data_dp) else None
    mkdir(sources_dp) if not isdir(sources_dp) else None
    mkdir(preprocessed_dp) if not isdir(preprocessed_dp) else None
    mkdir(output_dp) if not isdir(output_dp) else None
    
    # load sources' urls
    sources_urls = RawConfigParser()
    sources_urls.read(sources_fp)

    # download uniprot source data

    download_uniprot_files(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download reactome source data

    download_reactome_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download cellosaurus source data
    download_cellosaurus_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download hpa source data
    download_hpa_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download ctd source data
    download_ctd_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download phosphosite source data
    download_phosphosite_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download intact source data
    download_intact_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download sider source data
    download_sider_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download sider source data
    download_medgen_data(sources_dp=sources_dp, srcs_cp=sources_urls)
    # download drugbank source data
    if len(sys.argv) >= 3:
        db_user = sys.argv[1]
        db_pass = sys.argv[2]
        download_drugbank_data(sources_dp=sources_dp, srcs_cp=sources_urls, username=db_user, password=db_pass)
    elif 'DB_USER' in environ and 'DB_PASS' in environ:
        db_user = environ['DB_USER']
        db_pass = environ['DB_PASS']
        download_drugbank_data(sources_dp=sources_dp, srcs_cp=sources_urls, username=db_user, password=db_pass)
    # download kegg source data
    download_kegg_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download mesh source data
    download_mesh_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download Cutillas 20 data
    download_cutillas20_data(sources_dp=sources_dp, srcs_cp=sources_urls)

    # download SMPDB data
    download_smpdb_data(sources_dp=sources_dp, srcs_cp=sources_urls)
    # ----------------------------------------------------------------------
    # processing uniprot entries file
    uniprot_parser = UniProtTxtParser()
    uniprot_dp = join(preprocessed_dp, 'uniprot')
    mkdir(uniprot_dp) if not isdir(uniprot_dp) else None
    uniprot_output_fps = [join(uniprot_dp, fn) for fn in uniprot_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in uniprot_output_fps]))

    if invalid_md5:
        uniprot_parser.parse(sources_dp, uniprot_dp)
        for ofp in uniprot_output_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Uniprot processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing HPA entries file
    hpa_parser = HumanProteinAtlasParser()
    hpa_dp = join(preprocessed_dp, 'hpa')
    mkdir(hpa_dp) if not isdir(hpa_dp) else None
    hpa_files = ["hpa_antibodies.txt", "hpa_cellines_exp.txt", "hpa_tissues_exp.txt"]
    hpa_fps = [join(hpa_dp, fn) for fn in hpa_files]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hpa_fps]))
    if invalid_md5:
        hpa_parser.parse_database_xml(join(sources_dp, "proteinatlas.xml.gz"), hpa_dp)
        for ofp in hpa_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "HPA processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    # ----------------------------------------------------------------------
    # processing Cellosaurus database file
    cell_parser = CellosaurusParser()
    cell_dp = join(preprocessed_dp, 'cellosaurus')
    mkdir(cell_dp) if not isdir(cell_dp) else None
    cell_parser.parse_db_file(join(sources_dp, "cellosaurus_data.txt"), cell_dp)

    # TODO: export md5 hashes of resulting files as in other databases

    # ----------------------------------------------------------------------
    # processing DrugBank entries file
    drugbank_source_fp = join(sources_dp, "drugbank_all_full_database.xml.zip")
    drugbank_parser = DrugBankParser()
    db_dp = join(preprocessed_dp, 'drugbank')
    mkdir(db_dp) if not isdir(db_dp) else None
    drugbank_fps = [join(db_dp, fn) for fn in drugbank_parser.filelist]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in drugbank_fps]))
    if invalid_md5:
        # Skip drugbank processing if the file is not in the source folder
        if exists(drugbank_source_fp):
            drugbank_parser.parse_drugbank_xml(drugbank_source_fp, db_dp)
            for ofp in drugbank_fps:
                export_file_md5(ofp)
        else:
            print(fail_sym + "Drugbank source not available >>> Skipping Drugbank processing")
    else:
        print(inf_sym + "DrugBank processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing MESH
    mesh_parser = MESHParser()
    mesh_diseases_fp = join(sources_dp, 'mesh_diseases.txt')
    mesh_supp_fp = join(sources_dp, "mesh_supp_concepts.xml")
    mesh_dp = join(preprocessed_dp, 'mesh')
    mkdir(mesh_dp) if not isdir(mesh_dp) else None
    mesh_fps = [join(mesh_dp, fn) for fn in mesh_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in mesh_fps]))
    if invalid_md5:
        mesh_parser.parse_mesh(mesh_diseases_fp, mesh_supp_fp, mesh_dp)
        for ofp in mesh_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "MESH processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    # ----------------------------------------------------------------------
    # processing KEGG links
    kegg_parser = KeggParser()
    kegg_diseases_fp = join(sources_dp, 'diseases.txt')
    kegg_dp = join(preprocessed_dp, 'kegg')
    mkdir(kegg_dp) if not isdir(kegg_dp) else None
    kegg_fps = [join(kegg_dp, fn) for fn in kegg_parser.filelist]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in kegg_fps]))
    if invalid_md5:
        kegg_parser.parse_kegg(kegg_diseases_fp, kegg_dp)
        for ofp in kegg_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "KEGG processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing Reactome entries file
    reactome_parser = ReactomeParser()
    reactome_dp = join(preprocessed_dp, 'reactome')
    mkdir(reactome_dp) if not isdir(reactome_dp) else None
    reactome_fps = [join(reactome_dp, fn) for fn in reactome_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in reactome_fps]))

    if invalid_md5:
        reactome_parser.parse_reactome(sources_dp, reactome_dp)
        for ofp in reactome_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Reactome processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing CTD entries file
    ctd_parser = CTDParser()
    ctd_dp = join(preprocessed_dp, 'ctd')
    mkdir(ctd_dp) if not isdir(ctd_dp) else None
    ctd_fps = [join(ctd_dp, fn) for fn in ctd_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in ctd_fps]))

    if invalid_md5:
        ctd_parser.parse_ctd(sources_dp, join(uniprot_dp, 'uniprot_metadata.txt'), ctd_dp)
        for ofp in ctd_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "CTD processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    
    # ----------------------------------------------------------------------
    # processing Phosphosite entries file
    phosphosite_parser = PhosphositeParser()
    phosphosite_dp = join(preprocessed_dp, 'phosphosite')
    mkdir(phosphosite_dp) if not isdir(phosphosite_dp) else None
    phosphosite_fps = [join(phosphosite_dp, fn) for fn in phosphosite_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in phosphosite_fps]))

    if invalid_md5:
        phosphosite_parser.parse_phosphosite(sources_dp, phosphosite_dp)
        for ofp in phosphosite_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "PhosphoSitePlus processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing Intact zip file
    intact_parser = IntactParser()
    intact_dp = join(preprocessed_dp, 'intact')
    mkdir(intact_dp) if not isdir(intact_dp) else None
    intact_fps = [join(intact_dp, fn) for fn in intact_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in intact_fps]))

    if invalid_md5:
        intact_parser.parse_intact(sources_dp, intact_dp, join(uniprot_dp, 'uniprot_ppi.txt'))
        for ofp in intact_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Intact processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing Sider files
    sider_parser = SiderParser()
    sider_dp = join(preprocessed_dp, 'sider')
    mkdir(sider_dp) if not isdir(sider_dp) else None
    sider_fps = [join(sider_dp, fn) for fn in sider_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in sider_fps]))

    if invalid_md5:
        sider_parser.parse_sider(sources_dp, sider_dp)
        for ofp in sider_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Sider processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    
    # ----------------------------------------------------------------------
    # processing MedGen files
    medgen_parser = MedgenParser()
    medgen_dp = join(preprocessed_dp, 'medgen')
    mkdir(medgen_dp) if not isdir(medgen_dp) else None
    medgen_fps = [join(medgen_dp, fn) for fn in medgen_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in medgen_fps]))

    if invalid_md5:
        medgen_parser.parse_medgen(sources_dp, medgen_dp)
        for ofp in medgen_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "MedGen processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing Cutillas20 files
    cutillas20_parser = Cutillas20Parser()
    cutillas20_dp = join(preprocessed_dp, 'cutillas20')
    mkdir(cutillas20_dp) if not isdir(cutillas20_dp) else None
    cutillas20_fps = [join(cutillas20_dp, fn) for fn in cutillas20_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in cutillas20_fps]))

    if invalid_md5:
        cutillas20_parser.parse_phosphorylation(sources_dp, cutillas20_dp)
        for ofp in cutillas20_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Cutillas20 processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing SMPDB files
    smpdb_parser = SmpdbParser()
    smpdb_dp = join(preprocessed_dp, 'smpdb')
    mkdir(smpdb_dp) if not isdir(smpdb_dp) else None
    smpdb_fps = [join(smpdb_dp, fn) for fn in smpdb_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in smpdb_fps]))

    if invalid_md5:
        smpdb_parser.parse_pathways(sources_dp, smpdb_dp)
        for ofp in smpdb_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "SMPDB processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)


if __name__ == '__main__':
    preprocess_graph()
