# -*- coding: utf-8 -*-
import gzip
import re
import requests
from os.path import join
from ..util.extras import *
from timeit import default_timer as timer
import xml.etree.ElementTree as ET
from zipfile import ZipFile


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


def sanatize_text(text):
    """ Replace non alphanumeric characters in text with '_'

    Parameters
    ----------
    text : str
        text to sanatize

    Returns
    -------
    text
        the sanatized text
    """
    if text is None:
        return text
    return re.sub('[^a-zA-Z0-9]', '_', text.strip())


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



class SetWriter():
    """
    Utility class for writing DrugBank statements
    Enforces uniqueness of statements between written between flushes
    Set clear_on_flush to false to enforce uniquness for on all writes
    (should not be set for very large datasets)
    """
    def __init__(self, path):
        """
        Initialize a new SetWriter

        Parameters
        ----------
        """
        self._lines = []
        self._lineset = set()
        self._fd = open(path, 'w')
        self._clear_on_flush = True
        self._closed = False


    @property
    def clear_on_flush(self):
        return self._clear_on_flush


    @clear_on_flush.setter
    def clear_on_flush(self, value):
        self._clear_on_flush = value 


    def write(self, line):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        
        if line in self._lineset:
            return
        self._lineset.add(line)
        self._lines.append(line)


    def flush(self):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        self._fd.writelines(self._lines)
        self._lines = []
        if self._clear_on_flush:
            self._lineset = set()


    def close(self):
        if len(self._lines) > 0:
            self.flush()
        self._lineset = set()
        self._fd.close()


class DrugBankParser():
    """
    A DrugBank data parser
    """
    def __init__(self):
        """
        """
        self._filemap = {
            "interaction" : "db_interactions.txt",
            "target" : "db_targets.txt",
            "pathway" : "db_pathways.txt",
            "meta" : "db_meta.txt",
            "mesh" : "db_mesh.txt",
            "classification" : "db_classification.txt",
            "atc" : "db_atc.txt"
        }
        self._ns = {'db':'http://www.drugbank.ca'}

    @property
    def filelist(self):
        return [s for s in self._filemap.values()]

    def __parse_target(self, target_element, drug_id, rel_type, output_fd):
        """
        Parse a drug target
        
            targets with actions are set as unknown
            targets with action 'other' are set as unknown
            targets with action 'other/unknown' are set as unknown

        Parameters
        ----------
        target_element : xml.etree.ElementTree.Element
            xml element
        drug_id : string
            id of the drug
        rel_type : string
            type of target
        output : SetWriter
            writer for statements
        """

        #Link to uniprot
        poly = target_element.find('./db:polypeptide', self._ns)
        if poly is None:
            return
        poly_id = None
        if 'source' in poly.attrib and poly.attrib['source'] == 'Swiss-Prot':
            poly_id  = poly.attrib['id']
        else:
            for extern_id in poly.findall('./db:external-identifiers/db:external-identifier', self._ns):
                res = extern_id.find('db:resource', self._ns)
                val = extern_id.find('db:identifier', self._ns)
                if res is not None and val is not None and res.text=='UniprotKB':
                    poly_id = sanatize_text(val.text)
                    break
        if poly_id is None:
            return

        #gather any references to pubmed
        pmids = target_element.findall('./db:references/db:articles/db:article/db:pubmed-id', self._ns)
        pmids = [sanatize_text(pmid.text) for pmid in filter(lambda x: x.text is not None, pmids)]
        ref_string = ''
        if len(pmids) > 0:
            ref_string = ','.join(pmids)

        #gather all actions
        actions = target_element.findall('./db:actions/db:action', self._ns)
        formatted_actions = []
        for action in actions:
            action_text = sanatize_text(action.text)
            #
            if action_text == 'other' or action_text == 'other_unknown':
                action_text = 'unknown'
            if action_text == '':
                continue
            formatted_actions.append(action_text)
       
        #If no action provided set it to unknown
        if len(formatted_actions) == 0:
            formatted_actions = ['unknown']


        #create an extended quad for each action including references
        for action in formatted_actions:
            if len(pmids) > 0:
                output_fd.write(f'{drug_id}\t{rel_type}\t{poly_id}\t{action}\t{ref_string}\n')
            else:
                output_fd.write(f'{drug_id}\t{rel_type}\t{poly_id}\t{action}\n')


    def __parse_drug_interaction(self, interaction_element, drug_id, output):
        """
        Parse a drug interaction

        Parameters
        ----------
        interaction_element : xml.etree.ElementTree.Element
            xml element
        drug_id : string
            id of the drug
        output : SetWriter
            writer for statements
        """
        dest = interaction_element.find('./db:drugbank-id', self._ns)
        if dest is None:
            raise Exception('Interaction does not contain destination')
        dest_text = sanatize_text(dest.text)
        if dest_text is not None and dest_text != '':
            output.write(f'{drug_id}\tDRUG_INTERACTION\t{dest_text}\n')


    def __parse_atc_code(self, code_element, drug_id, output):
        """
        Parse a drug atc codes 
        the code string encodes all levels of the code hierarchy
        for example B01AE02 
        B       : C1 
        B01     : C2
        B01A    : C3
        B01AE   : C4
        B01AE02 : C5

        Parameters
        ----------
        code_element : xml.etree.ElementTree.Element
            xml element
        drug_id : string
            id of the drug
        output : SetWriter
            writer for statements
        """
        code = code_element.get('code')
        output.write(f'{drug_id}\tDRUG_ATC_C1\tATC:{code[0:1]}\n')
        output.write(f'{drug_id}\tDRUG_ATC_C2\tATC:{code[0:3]}\n')
        output.write(f'{drug_id}\tDRUG_ATC_C3\tATC:{code[0:4]}\n')
        output.write(f'{drug_id}\tDRUG_ATC_C4\tATC:{code[0:5]}\n')
        output.write(f'{drug_id}\tDRUG_ATC_C5\tATC:{code}\n')


    def __parse_pathway(self, pathway_element, drug_id, output):
        """
        Parse a drug pathway

        Parameters
        ----------
        pathway_element : xml.etree.ElementTree.Element
            xml element
        drug_id : string
            id of the drug
        output : SetWriter
            writer for statements
        """
        pid = pathway_element.find('./db:smpdb-id', self._ns)
        if pid is None:
            return
        pid = pid.text
        output.write(f'{drug_id}\tDRUG_PATHWAY\t{pid}\n')
        category = pathway_element.find('./db:category', self._ns)
        if category is not None:
            category_text = sanatize_text(category.text)
            if category_text is not None and category_text != '':
                output.write(f'{pid}\tPATHWAY_CATEGORY\t{category.text.strip()}\n')
        for enzyme in pathway_element.findall('./db:enzymes/db:uniprot-id', self._ns):
            enzyme_text = sanatize_text(enzyme.text)
            if enzyme_text is not None and enzyme_text != '':
                output.write(f'{pid}\tPATHWAY_ENZYME\t{enzyme_text}\n')


    def __parse_drug(self, drug_element, output_writers):
        """
        Parse a top level xml drug entry
        
        Parameters
        ----------
        drug_element : xml.etree.ElementTree.Element
            xml element
        output_writers: dict
            maps section names to their SetWriters
        """

        #
        #Parse drug metadata
        meta_fd = output_writers['meta']
        drug_id_elem = drug_element.find('./db:drugbank-id[@primary="true"]', self._ns)
        
        if drug_id_elem is None:
            raise Exception('Primary id not found')
        
        drug_id = drug_id_elem.text
        meta_fd.write(f'{drug_id}\tTYPE\tDRUG\n')
        name = drug_element.find('./db:name', self._ns)
        if name is not None:
            name_text = sanatize_text(name.text)
            if name_text is not None and name_text != '':
                meta_fd.write(f'{drug_id}\tNAME\t{name_text}\n')

        for synonym in drug_element.findall('./db:synonyms/db:synonym[@language="english"]', self._ns):
            syn_text = sanatize_text(synonym.text)
            if syn_text is not None and syn_text != '':
                meta_fd.write(f'{drug_id}\tSYNONYM\t{syn_text}\n')

        for group in drug_element.findall('./db:groups/db:group', self._ns):
            group_text = sanatize_text(group.text)
            if group_text is not None and group_text != '':
                meta_fd.write(f'{drug_id}\tGROUP\t{group_text}\n')

        for pmid in drug_element.findall('./db:general-references/db:articles/db:article/db:pubmed-id', self._ns):
            pmid_text = sanatize_text(pmid.text)
            if pmid_text is not None and pmid_text != '':
                meta_fd.write(f'{drug_id}\tPUBMED_ARTICLE\t{pmid_text}\n')

        for product in drug_element.findall('./db:products/db:product/db:name', self._ns):        
            product_text = sanatize_text(product.text)
            if product_text is not None and product_text != '':
                meta_fd.write(f'{drug_id}\tPRODUCT\t{product_text}\n')

        #
        #Parse drug classification
        classification_fd = output_writers['classification']
        classification = drug_element.find('./db:classification', self._ns)
        if classification is not None:
            for child in classification:
                if child.tag == '{%s}description' % self._ns['db']:
                    continue
                if child.text is not None and child.text != '':
                    c_type = child.tag.split('}')[-1]
                    value = sanatize_text(child.text)
                    if value is not None and value != '':
                        classification_fd.write(f'{drug_id}\t{c_type}\t{value}\n')

        #
        #Parse drug targets
        target_fd = output_writers['target']
        for target in drug_element.findall('./db:targets/db:target', self._ns):
            self.__parse_target(target, drug_id, 'DRUG_TARGET', target_fd)

        for carrier in drug_element.findall('./db:carriers/db:carrier', self._ns):
            self.__parse_target(carrier, drug_id, 'DRUG_CARRIER', target_fd)

        for transporter in drug_element.findall('./db:transporters/db:transporter', self._ns):
            self.__parse_target(transporter, drug_id, 'DRUG_TRANSPORTER', target_fd)

        for enzyme in drug_element.findall('./db:enzymes/db:enzyme', self._ns):
            self.__parse_target(enzyme, drug_id, 'DRUG_ENZYME', target_fd)

        #
        #Parse drug interactions
        interaction_fd = output_writers['interaction']
        for interaction in drug_element.findall('./db:drug-interactions/db:drug-interaction', self._ns):
            self.__parse_drug_interaction(interaction, drug_id, interaction_fd)

        #
        #Parse drug atc code categories
        atc_fd = output_writers['atc']
        for atc_code in drug_element.findall('./db:atc-codes/db:atc-code', self._ns):
            self.__parse_atc_code(atc_code, drug_id, atc_fd)

        #
        #Parse mesh categories
        mesh_fd = output_writers['mesh']
        for mesh_id in drug_element.findall('./db:categories/db:category/db:mesh-id', self._ns):
            mesh_id_text = sanatize_text(mesh_id.text)
            if mesh_id_text is not None and mesh_id_text != '':
                mesh_fd.write(f'{drug_id}\tMESH_CATEGORY\t{mesh_id.text}\n')

        #
        #Parse drug pathways
        pathway_fd = output_writers['pathway']
        for pathway in drug_element.findall('./db:pathways/db:pathway', self._ns):
            self.__parse_pathway(pathway, drug_id, pathway_fd)

        
    def parse_drugbank_xml(self, filepath, output_dp, filename='full database.xml'):
        """ Parse Drugbank xml file

        Parameters
        ----------
        filepath : str
            absolute file path of the drugbank zip file
        output_dp : str
            path of the output directory
        filename : str
            name of the xml file in the drugbank zip (default "full database.xml")
        """
        output_writers = {key: SetWriter(join(output_dp, fn)) for key, fn in self._filemap.items()}
        output_writers['pathway'].clear_on_flush = False

        with ZipFile(filepath, 'r') as dbzip:
            with dbzip.open(filename, force_zip64=True) as xmlfile:
                print_section_header("Parsing Drugbank XML file (%s)" % 
                    (bcolors.OKGREEN + filepath + "/" + filename + bcolors.ENDC))
                start = timer()
                nb_entries = 0
                for event, elem in ET.iterparse(xmlfile):
                    #Check the length of the drug element as pathways also contain drug elements
                    if elem.tag=='{http://www.drugbank.ca}drug' and len(elem) > 2:
                        nb_entries += 1
                        if nb_entries % 5 == 0:
                            speed = nb_entries / (timer() - start)
                            msg = prc_sym + "Processed (%d) entries.  Speed: (%1.5f) entries/second" % (nb_entries, speed)
                            print("\r" + msg, end="", flush=True)
                        self.__parse_drug(elem, output_writers)
                        elem.clear()

                        #Flush the output buffers
                        for writer in output_writers.values():
                            writer.flush()
                print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)
        
        for writer in output_writers.values():
            writer.close()
