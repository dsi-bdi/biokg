# -*- coding: utf-8 -*-
import gzip
import re
import requests
from os.path import join
from ..util.extras import *
from timeit import default_timer as timer
import xml.etree.ElementTree as ET


P_UNIPROT_CODE = re.compile("[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]")
P_DISEASE_CODE = re.compile("MIM:\d+")
P_GO_ANOT_CODE = re.compile("GO:\d{7}")
P_PFAM____CODE = re.compile("PF:\d{3,6}")
P_MODRES__CODE = re.compile("^\w*(\.|\;)")
P_EC_NUM__CODE = re.compile("EC=\d+.\d+.\d+.\d+")
PUBMED_ID_CODE = re.compile("PubMed=\d+")
SEQ_RANGE_CODE = re.compile("\w+\s+\d+\.\.\d+")
SEQ_NOTE__CODE = re.compile("\w+\s+\d+\s")


def export_quads(triplets, file_descriptor):
    """ Export quads to a file

    Parameters
    ----------
    triplets : list
        list of quads
    file_descriptor : file
        file descriptor of an open writable file
    """

    for s, p, o, q in triplets:
        file_descriptor.write("%s\t%s\t%s\t%s\n" % (s, p, o, q))


def export_triplets(triplets, file_descriptor):
    """ Export triplets to a file
    
    Parameters
    ----------
    triplets : list
        list of triplets
    file_descriptor : file
        file descriptor of an open writable file

    """

    for s, p, o in triplets:
        file_descriptor.write("%s\t%s\t%s\n" % (s, p, o))


class UniProtTxtParser:
    """
    A UNIPROT database text file parser class
    """
    def __init__(self):
        """

        """
        self.filepath = ""
        self.output_dp = ""

    def __parse_txt_entry(self, entry):
        """ Process a Uniprot txt entry

        Parameters
        ----------
        entry : list
            list of str lines representing the entry lines

        Returns
        -------
        list
            list of extracted facts
        """
        entry_facts = []
        entry_metadata = []
        entry_ppi = []

        entry_dictionary = dict()
        for line in entry:
            line_prefix = line[:2]
            if line_prefix == "  ":
                line_prefix = "AA"

            if line_prefix not in entry_dictionary:
                entry_dictionary[line_prefix] = ""

            entry_dictionary[line_prefix] += " " + line[2:].lstrip()

        # ------------------------------------------------------------------------
        # Processing IDs and UniProt ACs
        # ------------------------------------------------------------------------
        entry_code = entry_dictionary["AC"].strip().split(";")[0]
        row = entry_dictionary["AC"].split()
        if len(row) > 2:
            for idx in range(2, len(row)):
                entry_metadata.append([entry_code, "OTHER_ID", row[idx][:-1]])

        entry_name = entry_dictionary["ID"].strip().split(" ")[0].split("_")[0]
        entry_species = entry_dictionary["ID"].strip().split(" ")[0].split("_")[-1]

        entry_metadata.append([entry_code, "NAME", entry_name])
        entry_metadata.append([entry_code, "SPECIES", entry_species])

        # ------------------------------------------------------------------------
        # Processing DE prefix section
        # ------------------------------------------------------------------------
        if "DE" in entry_dictionary:
            full_names = [v[:-1].strip().replace("Full=", "") for v in
                          re.findall("Full=[^;]+;", entry_dictionary["DE"])]
            short_names = [v[:-1].strip().replace("Short=", "") for v in
                           re.findall("Short=[^;]+;", entry_dictionary["DE"])]

            for full_n in full_names:
                entry_metadata.append([entry_code, "FULL_NAME", full_n])
            for short_n in short_names:
                entry_metadata.append([entry_code, "SHORT_NAME", short_n])

        # ------------------------------------------------------------------------
        # Processing OC prefix section
        # ------------------------------------------------------------------------
        if "OC" in entry_dictionary:
            organism_classes = [c.strip() for c in entry_dictionary["OC"].strip().split(";")]
            for oc in organism_classes:
                entry_metadata.append([entry_code, "ORGANISM_CLASS", oc])

        # ------------------------------------------------------------------------
        # Processing RX prefix section
        # ------------------------------------------------------------------------
        if "RX" in entry_dictionary:
            pubmed_ids = ["pubmed:"+v.split("=")[-1] for v in re.findall(PUBMED_ID_CODE, entry_dictionary["RX"])]
            for pm in pubmed_ids:
                entry_metadata.append([entry_code, "RELATED_PUBMED_ID", pm])

        # ------------------------------------------------------------------------
        # Processing KW prefix section
        # ------------------------------------------------------------------------
        if "KW" in entry_dictionary:
            keywords = [v.strip() for v in entry_dictionary["KW"].replace(".", "").strip().split(";")]
            for kw in keywords:
                entry_metadata.append([entry_code, "RELATED_KEYWORD", kw])

        # ------------------------------------------------------------------------
        # Processing DR prefix section
        # ------------------------------------------------------------------------
        if "DR" in entry_dictionary:
            links_lines = entry_dictionary["DR"].strip().split(".")
            links_lines_dict = dict()
            for line in links_lines:
                db_name = line.strip().split(";")[0]
                if db_name not in links_lines_dict:
                    links_lines_dict[db_name] = []
                links_lines_dict[db_name].append(line.strip())

            if "GO" in links_lines_dict:
                go_lines = links_lines_dict["GO"]
                for line in go_lines:
                    go_code = line.split(";")[1].strip()
                    if "; F:" in line:
                        go_code_type = "GO_BIO_FUNC"
                    elif "; P:" in line:
                        go_code_type = "GO_MOL_PROC"
                    else:
                        go_code_type = "GO_CELL_LOC"
                    entry_facts.append([entry_code, go_code_type, go_code])

            if "HPA" in links_lines_dict:
                hpa_lines = links_lines_dict["HPA"]
                for line in hpa_lines:
                    hpa_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "RELATED_ANTIBODY", hpa_code])

            if "Reactome" in links_lines_dict:
                reactome_lines = links_lines_dict["Reactome"]
                for line in reactome_lines:
                    reactome_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "RELATED_PATHWAY", reactome_code])

            if "DrugBank" in links_lines_dict:
                drugbank_lines = links_lines_dict["DrugBank"]
                for line in drugbank_lines:
                    drugbank_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "TARGET_OF_DRUG", drugbank_code])

            if "InterPro" in links_lines_dict:
                interpro_lines = links_lines_dict["InterPro"]
                for line in interpro_lines:
                    interpro_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "SEQ_ANNOTATION", interpro_code])

        # ------------------------------------------------------------------------
        # Processing OC prefix section
        # ------------------------------------------------------------------------
        if "CC" in entry_dictionary:
            comments_cats_dict = dict()
            comments_list = [v.strip() for v in entry_dictionary["CC"].strip().split("-!-") if len(v)>3]
            for comment in comments_list:
                comment_cat = comment[:comment.find(":")].strip()
                comment_val = comment[comment.find(":")+1:].strip()
                comments_cats_dict[comment_cat] = comment_val
            if "INTERACTION" in comments_cats_dict:
                interactors_uniprot_acs = re.findall(P_UNIPROT_CODE, comments_cats_dict["INTERACTION"])
                for up_id in interactors_uniprot_acs:
                    entry_ppi.append([entry_code, "INTERACTS_WITH", up_id])
            if "DISEASE" in comments_cats_dict:
                disease_codes = re.findall(P_DISEASE_CODE, comments_cats_dict["DISEASE"])
                for c in disease_codes:
                    entry_facts.append([entry_code, "RELATED_DISEASE", c])

        # ------------------------------------------------------------------------
        # Processing FT prefix section [Sequence annotations]
        # ------------------------------------------------------------------------
        if "FT" in entry_dictionary:
            ft_content = entry_dictionary["FT"]
            seq_ranges = [v.strip().split() for v in re.findall(SEQ_RANGE_CODE, ft_content)]
            seq_annotations = [v.strip().split() for v in re.findall(SEQ_NOTE__CODE, ft_content)]

        # ------------------------------------------------------------------------
        return entry_facts, entry_metadata, entry_ppi

    def parse(self, filepath, output_dp):
        """ Parse a Uniprot textual data file and output findings to a set of files in a specified directory

        Parameters
        ----------
        filepath : str
            absolute path to a Uniprot textual data file

        output_dp : str
            absolute path of the output directory
        """
        self.filepath = filepath
        self.output_dp = output_dp

        facts_fd = open(join(output_dp, "uniprot_facts.txt"), "w")
        metadata_fd = open(join(output_dp, "uniprot_metadata.txt"), "w")
        ppi_fd = open(join(output_dp, "uniprot_ppi.txt"), "w")

        line_index = 0
        nb_facts = 0
        nb_metadata = 0
        nb_ppi = 0
        with gzip.open(filepath, 'rt') as fd:
            current_entry = []
            print_section_header("Parsing Uniprot file (%s)" % (bcolors.OKGREEN + filepath + bcolors.ENDC))
            start = timer()
            eof = False
            while not eof:
                raw_line = fd.readline()
                line = raw_line.rstrip()
                if line != "//":
                    current_entry.append(line)
                else:
                    facts, metadata, ppi = self.__parse_txt_entry(current_entry)

                    nb_facts += len(facts)
                    nb_metadata += len(metadata)
                    nb_ppi += len(ppi)

                    export_triplets(facts, facts_fd)
                    export_triplets(metadata, metadata_fd)
                    export_triplets(ppi, ppi_fd)
                    facts_fd.flush()
                    metadata_fd.flush()
                    ppi_fd.flush()
                    current_entry = []

                line_index += 1
                if line_index % 5000 == 0:
                    speed = int(line_index / (timer()-start))
                    msg = prc_sym + "Processing (%d) lines/second => [facts:%d - metadata:%d - ppi:%d - total: %d]" \
                          % (speed, nb_facts, nb_metadata, nb_ppi, nb_facts + nb_metadata + nb_ppi)
                    print("\r" + msg, end="", flush=True)

                if raw_line == "":
                    eof = True
                    print(done_sym + " Took %1.2f Seconds." % (timer()-start), flush=True)

        facts_fd.close()
        metadata_fd.close()
        ppi_fd.close()


class HumanProteinAtlasParser:
    """
    A Human Protein Atlas database parser
    """
    def __init__(self):
        """
        initialise new class instance
        """
        self.filepath = ""
        self.output_dp = ""

    def __parse_hpa_xml_entry(self, entry):
        """ parse an xml element of an HPA entry

        Parameters
        ----------
        entry : xml.etree.ElementTree.Element
            xml element

        Returns
        -------
        dict
            dictionary of parsed data
        """
        parsed_data = {"rna": [], "tissue": [], "ab": []}

        antibody_list = entry.findall("antibody")
        tissue_expression_list = entry.findall("tissueExpression")
        rna_expression_list = entry.findall("rnaExpression")

        if tissue_expression_list is not None:
            for tissue_expression in tissue_expression_list:
                data_elements = tissue_expression.findall("data")
                for el in data_elements:
                    cell_line = el.find("cellLine")
                    if cell_line is not None:
                        cell_line = cell_line.text
                        cell_line_level = el.find("level").text
                        parsed_data["tissue"].append(["CELLINE", cell_line, cell_line_level])

                    tissue = el.find("tissue").text.replace(" ", "_")
                    tissue_level = el.find("level").text
                    t_cells = el.findall("tissueCell")
                    parsed_data["tissue"].append(["TISSUE", tissue, tissue_level])
                    for tc in t_cells:
                        cell_type = tc.find("cellType").text.replace(" ", "_")
                        cell_type_level = tc.find("level").text
                        parsed_data["tissue"].append(["TISSUE", tissue + "__" + cell_type, cell_type_level])
        if rna_expression_list is not None:
            for rna_expression in rna_expression_list:
                data_elements = rna_expression.findall("data")
                for el in data_elements:
                    cell_line = el.find("cellLine")
                    if cell_line is not None:
                        cell_line_name = cell_line.text.replace(" ", "_")
                        cell_line_level = "NormRNA_EXP:" + el[1].attrib["expRNA"]
                        cellosaurus_id = cell_line.attrib["cellosaurusID"]
                        cellosaurus_id = "NA" if cellosaurus_id == "" else cellosaurus_id
                        organ = cell_line.attrib["organ"]
                        suffix = "%s#%s" % (organ, cellosaurus_id)
                        parsed_data["rna"].append(["CELLINE", cell_line_name, cell_line_level, suffix])

                    tissue = el.find("tissue")
                    if tissue is not None:
                        tissue = tissue.text.replace(" ", "_")
                        tissue_level = el.find("level").text
                        t_cells = el.findall("tissueCell")
                        parsed_data["rna"].append(["RNA", tissue, tissue_level])
                        for tc in t_cells:
                            cell_type = tc.find("cellType").text.replace(" ", "_")
                            cell_type_level = tc.find("level").text
                            parsed_data["rna"].append(["RNA", tissue + "__" + cell_type, cell_type_level])

        if antibody_list is not None:
            for antibody_elem in antibody_list:
                ab_id = antibody_elem.attrib["id"]
                antigen = antibody_elem.find("antigenSequence").text
                parsed_data["ab"].append(["ANTIBODY", ab_id, antigen])
        return parsed_data

    def parse_database_xml(self, filepath, output_dp):
        """ Parse HPA xml file

        Parameters
        ----------
        filepath : str
            absolute file path of the hpa xml file
        output_dp : str
            path of the output directory
        """
        self.filepath = filepath
        self.output_dp = output_dp
        xml_fd = gzip.open(filepath)
        cl_exp_fp = join(output_dp, "hpa_cellines_exp.txt")
        tissue_exp_fp = join(output_dp, "hpa_tissues_exp.txt")
        ab_data_fp = join(output_dp, "hpa_antibodies.txt")

        cl_exp_fd = open(cl_exp_fp, "w")
        tissue_exp_fd = open(tissue_exp_fp, "w")
        ab_data_fd = open(ab_data_fp, "w")

        print_section_header("Parsing HPA XML file (%s)" % (bcolors.OKGREEN + filepath + bcolors.ENDC))
        start = timer()
        nb_entries = 0
        for event, entry in ET.iterparse(xml_fd, events=('start', 'end')):
            if entry.tag == "entry" and event == "end" and len(entry.getchildren()) > 2:
                nb_entries += 1
                if nb_entries % 5 == 0:
                    speed = nb_entries / (timer() - start)
                    msg = prc_sym + "Processed (%d) entries.  Speed: (%1.5f) entries/second" % (nb_entries, speed)
                    print("\r" + msg, end="", flush=True)

                id_elem = entry.find("identifier")
                if len(id_elem) >= 1:
                    db_name = id_elem[0].attrib["db"]
                    if db_name == "Uniprot/SWISSPROT":
                        entry_id = id_elem[0].attrib["id"]
                        entry_data = self.__parse_hpa_xml_entry(entry)
                        # export antibody data
                        for _, ab, ag in entry_data["ab"]:
                            ag = ag if ag is not None else "-"
                            ab_data_fd.write("%s\t%s\t%s\n" % (ab, entry_id, ag))
                        # export tissue data
                        for context, ts, level in entry_data["tissue"]:
                            if level is not None:
                                if context == "TISSUE":
                                    tissue_exp_fd.write("%s\t%s\t%s\n" % (entry_id, ts, level))
                                else:
                                    cl_exp_fd.write("%s\t%s\t%s\n" % (entry_id, ts, level))
                        # export rna data
                        for rna_data in entry_data["rna"]:
                            if len(rna_data) == 3:
                                _, cl, level = rna_data
                                tissue_exp_fd.write("%s\t%s\t%s\n" % (entry_id, cl, level))
                            if len(rna_data) == 4:
                                _, cl, level, organ = rna_data
                                cl_exp_fd.write("%s\t%s\t%s\t%s\n" % (entry_id, cl, organ, level))

        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)
        ab_data_fd.close()
        tissue_exp_fd.close()
        cl_exp_fd.close()
