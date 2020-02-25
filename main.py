# -*- coding: utf-8 -*-

from os.path import join, isdir
from os import mkdir
from configparser import RawConfigParser

from biokg.loader import *
from biokg.util.extras import program_header
from biokg.processing.parsers import UniProtTxtParser, HumanProteinAtlasParser,\
    KeggParser, DrugBankParser, ReactomeParser
from biokg.util.io import export_file_md5, file_has_valid_md5


def main():
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

    # ----------------------------------------------------------------------
    # processing uniprot entries file
    uniprot_parser = UniProtTxtParser()
    uniprot_entries_fp = join(sources_dp, "swissprot_human_entries.txt.gz")
    uniprot_output_files = ["uniprot_facts.txt", "uniprot_metadata.txt", "uniprot_ppi.txt"]
    uniprot_output_fps = [join(output_dp, fn) for fn in uniprot_output_files]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in uniprot_output_fps]))

    if invalid_md5:
        uniprot_parser.parse(uniprot_entries_fp, output_dp)
        for ofp in uniprot_output_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Uniprot processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing HPA entries file
    hpa_parser = HumanProteinAtlasParser()
    hpa_files = ["hpa_antibodies.txt", "hpa_cellines_exp.txt", "hpa_tissues_exp.txt"]
    hpa_fps = [join(output_dp, fn) for fn in hpa_files]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hpa_fps]))
    if invalid_md5:
        hpa_parser.parse_database_xml(join(sources_dp, "proteinatlas.xml.gz"), preprocessed_dp)
        for ofp in hpa_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "HPA processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing DrugBank entries file
    drugbank_parser = DrugBankParser()
    drugbank_fps = [join(output_dp, fn) for fn in drugbank_parser.filelist]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in drugbank_fps]))
    if invalid_md5:
        drugbank_parser.parse_drugbank_xml(join(sources_dp, "drugbank_all_full_database.xml.zip"), output_dp)
        for ofp in drugbank_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "DrugBank processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing KEGG links
    kegg_parser = KeggParser()
    kegg_fp = join(output_dp, kegg_parser.filename)
    invalid_md5 = not file_has_valid_md5(kegg_fp)
    if invalid_md5:
        kegg_parser.parse_kegg(output_dp)
        export_file_md5(kegg_fp)
    else:
        print(inf_sym + "KEGG processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

    # ----------------------------------------------------------------------
    # processing Reactome entries file
    reactome_parser = ReactomeParser()
    reactome_fps = [join(preprocessed_dp, fn) for fn in reactome_parser.filenames]
    invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in reactome_fps]))

    if invalid_md5:
        reactome_parser.parse_reactome(sources_dp, preprocessed_dp)
        for ofp in reactome_fps:
            export_file_md5(ofp)
    else:
        print(inf_sym + "Reactome processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)


if __name__ == '__main__':
    main()
