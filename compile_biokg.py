from biodblinker import KEGGLinker, SiderLinker, MESHLinker
from os import makedirs, listdir, remove, walk
from os.path import join, isdir
from shutil import copy
import gzip

kegg_linker = KEGGLinker()
sider_linker = SiderLinker()
mesh_linker = MESHLinker()
data_root = 'data/preprocessed'
output_root = 'data/output'
core_root = 'data/biokg'
links_root = join(output_root, 'links')
meta_root = join(output_root, 'metadata')
properties_root = join(output_root, 'properties')
drug_properties_root = join(properties_root, 'drug')
protein_properties_root = join(properties_root, 'protein')
pathway_properties_root = join(properties_root, 'pathway')
disease_properties_root = join(properties_root, 'disease')
cell_properties_root = join(properties_root, 'cell')
mim_properties_root = join(properties_root, 'genetic_disorders')
other_root = join(output_root, 'other')

makedirs(output_root) if not isdir(output_root) else None
makedirs(links_root) if not isdir(links_root) else None
makedirs(meta_root) if not isdir(meta_root) else None
makedirs(properties_root) if not isdir(properties_root) else None
makedirs(drug_properties_root) if not isdir(drug_properties_root) else None
makedirs(protein_properties_root) if not isdir(protein_properties_root) else None
makedirs(pathway_properties_root) if not isdir(pathway_properties_root) else None
makedirs(disease_properties_root) if not isdir(disease_properties_root) else None
makedirs(cell_properties_root) if not isdir(cell_properties_root) else None
makedirs(mim_properties_root) if not isdir(mim_properties_root) else None
makedirs(other_root) if not isdir(other_root) else None
makedirs(core_root) if not isdir(core_root) else None

def get_all_proteins():
    """
    Get the set of uniprot proteins to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    uniprot_meta = join(data_root, 'uniprot', 'uniprot_metadata.txt')
    protein_set = set()
    with open(uniprot_meta, 'r') as fd:
        for line in fd:
            protein = line.split('\t')[0]
            protein_set.add(protein)

    return protein_set


def get_proteins_by_metadata(filter_field, accpeted_values):
    """
    Get the set of uniprot proteins for which the metadata value of field
    is ont of the accpeted_values

    Parameters
    ----------
    filter_field: str
        The metadata field to filter by
    accepted_values: list
        The set of accpetable values for the metadata field

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    uniprot_meta = join(data_root, 'uniprot', 'uniprot_metadata.txt')
    protein_set = set()
    with open(uniprot_meta, 'r') as fd:
        for line in fd:
            protein, field, value = line.strip().split('\t')
            if field == filter_field and value in accpeted_values:
                protein_set.add(protein)

    return protein_set


def get_all_drugs():
    """
    Get the set of drugs to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of drugs to include in the kg
    """
    db_meta = join(data_root, 'drugbank', 'db_meta.txt')

    db_drugs = set()

    with open(db_meta, 'r') as fd:
        for line in fd:
            drug = line.split('\t')[0]
            db_drugs.add(drug)

    return db_drugs


def get_all_mesh_diseases():
    """
    Get the set of mesh diseases to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of mesh diseases to include in the kg
    """
    files = [
        join(data_root, 'mesh', 'mesh_disease_meta.txt'),
        join(data_root, 'mesh', 'mesh_scr_disease_meta.txt')
    ]

    mesh_diseases = set()
    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                disease, meta, value = line.strip().split('\t')
                if meta == 'TYPE':
                    mesh_diseases.add(disease)

    return mesh_diseases


def get_all_genetic_disorders():
    """
    Get the set of uniprot proteins for which the metadata value of field
    is ont of the accpeted_values

    Parameters
    ----------
    filter_field: str
        The metadata field to filter by
    accepted_values: list
        The set of accpetable values for the metadata field

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    uniprot_meta = join(data_root, 'uniprot', 'uniprot_facts.txt')
    genetic_disorders = set()
    with open(uniprot_meta, 'r') as fd:
        for line in fd:
            protein, field, value = line.strip().split('\t')
            if field == 'RELATED_GENETIC_DISORDER':
                genetic_disorders.add(value)

    return genetic_disorders


def get_all_unique_ppi(protein_set):
    """
    Get the set of protein proteins interactions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of unique protein protein interactions
    """
    files = [
        join(data_root, 'uniprot', 'uniprot_ppi.txt'),
        join(data_root, 'reactome', 'reactome_ppi.txt'),
        join(data_root, 'intact', 'intact_ppi.txt')
    ]
    ppis = set()

    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                s, _, o = line.strip().split('\t')[:3]
                if s > o:
                    ppis.add((o, s))
                else:
                    ppis.add((s, o))
    unique_triples = []
    for s, o in ppis:
        if s in protein_set and o in protein_set:
            unique_triples.append((s, 'PPI', o))

    return unique_triples


def get_species_map():
    species_map = {}
    with open(join(data_root, 'uniprot', 'uniprot_metadata.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'SPECIES':
                species_map[s] = o

    return species_map


def write_ppi_by_species():
    """
    Get the set of protein proteins interactions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of unique protein protein interactions
    """
    species_map = get_species_map()
    unique_species = list(set(species_map.values()))
    species_root = join(links_root, 'PPI_SPECIES')
    makedirs(species_root) if not isdir(species_root) else None

    species_file_names = {}
    for species in set(species_map.values()):
        species_file_names[species] = join(species_root, f'{species}_ppi.txt')

    #species_files['OTHER'] = open(join(species_root, 'INTERSPECIES_ppi.txt'), 'w')
    
    files = [
        join(data_root, 'uniprot', 'uniprot_ppi.txt'),
        join(data_root, 'reactome', 'reactome_ppi.txt'),
        join(data_root, 'intact', 'intact_ppi.txt')
    ]
    ppis = set()

    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                s, _, o = line.strip().split('\t')[:3]
                if s > o:
                    ppis.add((o, s))
                else:
                    ppis.add((s, o))

    unique_triples = []
    species_files = {}
    for s, o in ppis:
        if s not in species_map or o not in species_map:
            continue
        s_species = species_map[s]
        o_species = species_map[o]
        if s_species == o_species:
            if s_species not in species_files:
                species_files[s_species] = open(species_file_names[s_species], 'w')
            species_files[s_species].write(f'{s}\tPPI\t{o}\n')
        else:
            if 'OTHER' not in species_files:
                species_files['OTHER'] = open(join(species_root, 'INTERSPECIES_ppi.txt'), 'w')
            species_files['OTHER'].write(f'{s}\tPPI\t{o}\n')

    for f in species_files.values():
        f.close()


def get_all_protein_sequence_annotations(protein_set):
    """
    Get the set of protein sequence annotations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of interpro protein sequence annotations
    list
        the list of prosite protein sequence annotations
    """
    seq_ann_root = join(protein_properties_root, 'sequence_annotations')
    makedirs(seq_ann_root) if not isdir(seq_ann_root) else None
    annotation_output_files = {
        'ACTIVE_SITE': open(join(seq_ann_root, 'protein_active_site.txt'), 'w'),
        'BINDING_SITE': open(join(seq_ann_root, 'protein_binding_site.txt'), 'w'),
        'CONSERVED_SITE': open(join(seq_ann_root, 'protein_conserved_site.txt'), 'w'),
        'DOMAIN': open(join(seq_ann_root, 'protein_domain.txt'), 'w'),
        'FAMILY': open(join(seq_ann_root, 'protein_family.txt'), 'w'),
        'HOMOLOGOUS_SUPERFAMILY': open(join(seq_ann_root, 'protein_homologous_superfamily.txt'), 'w'),
        'PTM': open(join(seq_ann_root, 'protein_ptm.txt'), 'w'),
        'REPEAT': open(join(seq_ann_root, 'protein_repeat.txt'), 'w'),
        'PS_SEQ_ANN': open(join(seq_ann_root, 'protein_prosite_sequence_annotation.txt'), 'w'),
        'GO_BP': open(join(protein_properties_root, 'protein_go_biological_process.txt'), 'w'),
        'GO_CC': open(join(protein_properties_root, 'protein_go_cellular_component.txt'), 'w'),
        'GO_MF': open(join(protein_properties_root, 'protein_go_molecular_function.txt'), 'w'),
        'RELATED_GENETIC_DISORDER': open(join(links_root, 'protein_genetic_disorders.txt'), 'w')
    }

    with open(join(data_root, 'uniprot', 'uniprot_facts.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p in annotation_output_files and s in protein_set:
                annotation_output_files[p].write(line)

    for fd in annotation_output_files.values():
        fd.close()


def get_all_protein_drug_interactions(protein_set, drug_set):
    """
    Get the set of protein drug interactions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein drug interactions
    drug_set: set
        the set of drugs used to filter the protein drug interactions

    Returns
    -------
    list
        the list of unique protein drug interactions
    """
    kegg_links = join(data_root, 'kegg', 'gene_drug.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    db_targets = join(data_root, 'drugbank', 'db_targets.txt')
    # ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    pdis = set()
    carriers = set()
    enzymes = set()
    transporters = set()
    targets = set()
    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([o])
            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(db_drugs[0]) == 0:
                continue
            if len(protein_acc[0]) == 0:
                continue
            for d in db_drugs[0]:
                for p in protein_acc[0]:
                    pdis.add((p, d))

    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'TARGET_OF_DRUG':
                pdis.add((s, o))

    with open(db_targets, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')[:3]
            pdis.add((o, s))

            if p == 'DRUG_CARRIER':
                carriers.add((o, s))
            elif p == 'DRUG_ENZYME':
                enzymes.add((o, s))
            elif p == 'DRUG_TRANSPORTER':
                transporters.add((o, s))
            else:
                targets.add((o, s))

    unique_triples = []
    function_triples = []
    for protein, drug in pdis:
        if protein in protein_set and drug in drug_set:
            unique_triples.append((drug, 'DPI', protein))

    for protein, drug in targets:
        if protein in protein_set and drug in drug_set:
            function_triples.append((drug, 'DRUG_TARGET', protein))

    for protein, drug in carriers:
        if protein in protein_set and drug in drug_set:
            function_triples.append((drug, 'DRUG_CARRIER', protein))

    for protein, drug in enzymes:
        if protein in protein_set and drug in drug_set:
            function_triples.append((drug, 'DRUG_ENZYME', protein))

    for protein, drug in transporters:
        if protein in protein_set and drug in drug_set:
            function_triples.append((drug, 'DRUG_TRANSPORTER', protein))

    return unique_triples, function_triples


def get_all_protein_disease_associations(protein_set, disease_set):
    """
    Get the set of protein disease associations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein disease associations
    disease_set: set
        the set of mesh diseases used to filter the protein disease associations

    Returns
    -------
    list
        the list of unique protein disease associations
    """
    kegg_links = join(data_root, 'kegg', 'gene_disease.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    ctd_protein_disease = join(data_root, 'ctd', 'ctd_protein_disease_association.txt')
    pdis = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            mesh_diseases = kegg_linker.convert_disease_to_mesh([o])
            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(protein_acc[0]) == 0:
                continue
            if len(mesh_diseases[0]) == 0:
                continue

            for p in protein_acc[0]:
                for d in mesh_diseases[0]:
                    pdis.add((p, d))

    with open(ctd_protein_disease, 'r') as fd:
        for line in fd:
            s, p, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                pdis.add((s, o))

    unique_triples = []
    for protein, disease in pdis:
        if protein in protein_set and disease in disease_set:
            unique_triples.append((protein, 'PROTEIN_DISEASE_ASSOCIATION', disease))

    return unique_triples


def get_all_protein_pathway_associations(protein_set):
    """
    Get the set of protein pathway associations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein pathway associations

    Returns
    -------
    list
        the list of unique protein pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'gene_pathway.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    ctd_pathways = [
        join(data_root, 'ctd', 'ctd_protein_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_protein_kegg_pathway_association.txt')
    ]

    db_pathways = join(data_root, 'drugbank', 'db_pathways.txt')

    ppis = set()
    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'RELATED_PATHWAY':
                ppis.add((s, o))

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')

            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(protein_acc[0]) == 0:
                protein_acc = [[s]]
            for acc in protein_acc[0]:
                ppis.add((acc, o))

    with open(db_pathways, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'PATHWAY_ENZYME':
                ppis.add((o, s))

    for f in ctd_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')
                if prov == 'CURATED':
                    ppis.add((s, o))

    unique_triples = []
    for protein, pathway in ppis:
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_PATHWAY_ASSOCIATION', pathway))

    return unique_triples


def get_all_drug_drug_interactions(drug_set):
    """
    Get the set of drug drug interactions to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug drug interactions

    Returns
    -------
    list
        the list of unique drug drug interactions
    """
    db_interactions = join(data_root, 'drugbank', 'db_ddi.txt')
    ddis = set()

    with open(db_interactions, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            if s > o:
                ddis.add((o, s))
            else:
                ddis.add((s, o))

    unique_triples = []
    for d1, d2 in ddis:
        if d1 in drug_set and d2 in drug_set:
            unique_triples.append((d1, 'DDI', d2))

    return unique_triples


def get_all_drug_pathway_associations(drug_set):
    """
    Get the set of drug pathway associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug pathway associations

    Returns
    -------
    list
        the list of unique drug pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'drug_pathway.txt')

    ctd_pathways = [
        join(data_root, 'ctd', 'ctd_drug_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_drug_kegg_pathway_association.txt')
    ]
    db_pathways = join(data_root, 'drugbank', 'db_pathways.txt')
    dpas = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([s])
            if len(db_drugs[0]) == 0:
                continue

            for d in db_drugs[0]:
                dpas.add((d, o))

    with open(db_pathways, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'DRUG_PATHWAY':
                dpas.add((s, o))

    for f in ctd_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')
                if prov == 'CURATED':
                    dpas.add((s, o))

    unique_triples = []
    for drug, pathway in dpas:
        if drug in drug_set:
            unique_triples.append((drug, 'DRUG_PATHWAY_ASSOCIATION', pathway))

    return unique_triples


def get_all_drug_disease_associations(drug_set, disease_set):
    """
    Get the set of drug pathway associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug pathway associations

    Returns
    -------
    list
        the list of unique drug pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'disease_drug.txt')
    ctd_links = join(data_root, 'ctd', 'ctd_drug_diseases.txt')
    ddis = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([o])
            mesh_diseases = kegg_linker.convert_disease_to_mesh([s])
            if len(db_drugs[0]) == 0:
                continue
            if len(mesh_diseases[0]) == 0:
                continue

            for dr in db_drugs[0]:
                for ds in mesh_diseases[0]:
                    ddis.add((dr, ds))

    with open(ctd_links, 'r') as fd:
        for line in fd:
            s, _, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                ddis.add((s, o))

    unique_triples = []
    for drug, disease in ddis:
        if drug in drug_set and disease in disease_set:
            unique_triples.append((drug, 'DRUG_DISEASE_ASSOCIATION', disease))
    return unique_triples


def get_all_drug_side_effects(drug_set):
    """
    Get the set of drug side effect associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug side effect associations

    Returns
    -------
    list
        the list of unique drug side effect associations
    """
    sider_effects = join(data_root, 'sider', 'sider_effects.txt')
    sider_indications = join(data_root, 'sider', 'sider_indications.txt')

    side_effects = set()
    indications = set()
    with open(sider_effects, 'r') as fd:
        for line in fd:
            d, _, o = line.strip().split('\t')
            db_ids = sider_linker.convert_drugs_to_drugbank([d])
            if len(db_ids[0]) == 0:
                continue

            for db_id in db_ids[0]:
                side_effects.add((db_id, o))

    with open(sider_indications, 'r') as fd:
        for line in fd:
            d, _, o = line.strip().split('\t')
            db_ids = sider_linker.convert_drugs_to_drugbank([d])
            if len(db_ids[0]) == 0:
                continue

            for db_id in db_ids[0]:
                indications.add((db_id, o))

    unique_triples = []
    for drug, effect in side_effects:
        if drug in drug_set:
            unique_triples.append((drug, 'DRUG_SIDEEFFECT_ASSOCIATION', effect))

    unique_indications = []
    for drug, indication in indications:
        if drug in drug_set:
            unique_indications.append((drug, 'DRUG_INDICATION_ASSOCIATION', indication))
    return unique_triples, unique_indications


def get_all_drug_atc_codes(drug_set):
    """
    Get the set of drug atc codes to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug drug interactions

    Returns
    -------
    list
        the list of unique drug atc codes
    """
    db_interactions = join(data_root, 'drugbank', 'db_atc.txt')
    datc = set()

    with open(db_interactions, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            datc.add((s, o))

    unique_triples = []
    for drug, atc in datc:
        if drug in drug_set:
            unique_triples.append((drug, 'DRUG_ATC_CODE', atc))

    return unique_triples


def get_all_disease_pathway_associations(disease_set):
    """
    Get the set of disease pathway associations to include in the kg

    Parameters
    ----------

    Returns
    -------
    list
        the list of unique disease pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'disease_pathway.txt')

    ddis = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')

            mesh_diseases = kegg_linker.convert_disease_to_mesh([s])
            if len(mesh_diseases[0]) == 0:
                continue

            for ds in mesh_diseases[0]:
                ddis.add((ds, o))

    unique_triples = []
    for disease, pathway in ddis:
        if disease in disease_set:
            unique_triples.append((disease, 'DISEASE_PATHWAY_ASSOCIATION', pathway))
    return unique_triples


def get_disease_tree(disease_set):
    """
    Get the set of disease pathway associations to include in the kg

    Parameters
    ----------

    Returns
    -------
    list
        the list of unique disease pathway associations
    """
    disease_tree = join(data_root, 'mesh', 'mesh_disease_tree.txt')
    distree = set()

    with open(disease_tree, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')           
            distree.add((s, o))

    unique_triples = []
    for disease, category in distree:
        if disease in disease_set:
            unique_triples.append((disease, 'DISEASE_SUPERGRP', category))
    return unique_triples


def get_pathway_rels():
    pathway_rels = set()
    with open(join(data_root, 'reactome', 'reactome_pathway_rels.txt'), 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')
            pathway_rels.add((s, o))

    unique_triples = []
    for pathway, parent in pathway_rels:
        unique_triples.append((pathway, 'HAS_PARENT_PATHWAY', parent))

    return unique_triples


def get_complex_pathway_rels():
    top_level_pathways = set()
    complex_pathways = set()
    with open(join(data_root, 'reactome', 'reactome_complex_pathway_rels.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'COMPLEX_PATHWAY':
                complex_pathways.add((s, o))
            else:
                top_level_pathways.add((s, o))

    unique_triples = []
    for _complex, pathway in complex_pathways:
        unique_triples.append((_complex, 'COMPLEX_IN_PATHWAY', pathway))

    tl_triples = []
    for _complex, pathway in top_level_pathways:
        tl_triples.append((_complex, 'COMPLEX_TOP_LEVEL_PATHWAY', pathway))
    return unique_triples, tl_triples


def get_protein_complex_rels(protein_set):
    protein_complex = set()
    with open(join(data_root, 'reactome', 'reactome_protein_complex_rels.txt'), 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            protein_complex.add((s, o))

    unique_triples = []
    for protein, _complex in protein_complex:
        if protein in protein_set:
            unique_triples.append((protein, 'MEMBER_OF_COMPLEX', _complex))

    return unique_triples


def get_all_protein_expressions(protein_set):
    """
    Get the set of protein tissue expressions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein tissue expressions
    Returns
    -------
    list
        the list of unique protein tissue expressions
    """
    hpa_tissue_expression = join(data_root, 'hpa', 'hpa_tissues_exp.txt')
    high_expression = set()
    medium_expression = set()
    low_expression = set()
    expression = set()
    tissue_structure = set()
    with open(hpa_tissue_expression, 'r') as fd:
        for line in fd:
            p, t, level = line.strip().split('\t')
            parts = t.split('__')
            if len(parts) == 2:
                tissue = parts[0]
                cell = parts[1]
                tissue_structure.add((cell, tissue))
            if level in ['low', 'medium', 'high']:
                expression.add((p, t))
            if level == 'low':
                low_expression.add((p, t))
            elif level == 'medium':
                medium_expression.add((p, t))
            elif level == 'high':
                high_expression.add((p, t))

    unique_triples = []
    # for protein, tissue in expression:
    #     if protein in protein_set:
    #         unique_triples.append((protein, 'Protein_Expression', tissue))
    hpa_other = join(other_root, 'hpa')
    makedirs(hpa_other) if not isdir(hpa_other) else None
    level_output = open(join(hpa_other, 'protein_expression_level.txt'), 'w')
    for protein, tissue in low_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tLOW\n')
    for protein, tissue in medium_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tMEDIUM\n')
    for protein, tissue in high_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tHIGH\n')

    structure_triples = []
    for cell, tissue in tissue_structure:
        structure_triples.append((cell, 'PART_OF_TISSUE', tissue))

    return unique_triples, structure_triples

def get_all_disease_genetic_disorder(disease_set, genetic_disorder_set):
    dgd = set()

    for disease in disease_set:
        mapped_mim = mesh_linker.convert_disease_to_omim([disease])
        for mim in mapped_mim[0]:
            title = f'MIM:{mim}'
            if title in genetic_disorder_set:
                dgd.add((disease, title))

    unique_triples = []
    for disease, disorder in dgd:
        unique_triples.append((disease, 'DISEASE_GENETIC_DISORDER', disorder))

    return unique_triples


def write_protein_cellline_expressions(protein_set):
    pcl = set()
    hpa_other = join(other_root, 'hpa')
    makedirs(hpa_other) if not isdir(hpa_other) else None
    with open(join(data_root, 'hpa', 'hpa_cellines_exp.txt'), 'r') as fd:
        for line in fd:
            pro, _, cl_tissue, exp = line.strip().split('\t')
            cl = cl_tissue.split('#')[1]
            if cl == 'NA':
                continue
            exp = exp.split(':')[1]
            pcl.add((pro, cl, exp))

    with gzip.open(join(hpa_other, 'protein_cellline_expression.txt.gz'), 'wt') as fd:
        for protein, cellline, expression in pcl:
            if protein in protein_set:
                fd.write(f'{protein}\t{cellline}\t{expression}\n')


def write_triples(triples, output_fp):
    with open(output_fp, 'w') as output:
        for s, p, o in triples:
            output.write(f'{s}\t{p}\t{o}\n')


def filter_ctd_drug_protein(protein_set):
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    ctd_drug_protein_filtered = join(links_root, 'ctd_drug_protein_interactions_SWISSPORT.txt' )
    with open(ctd_drug_protein_filtered, 'w') as output_fd:
        with open(ctd_drug_protein, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')[:4]
                if o in protein_set:
                    output_fd.write(line)


def write_uniprot_metadata():
    uniprot_meta_dp = join(meta_root, 'protein')
    makedirs(uniprot_meta_dp) if not isdir(uniprot_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(uniprot_meta_dp, 'uniprot_name.txt'), 'w'),
        'SHORT_NAME': open(join(uniprot_meta_dp, 'uniprot_shortname.txt'), 'w'),
        'FULL_NAME': open(join(uniprot_meta_dp, 'uniprot_fullname.txt'), 'w'),
        'ORGANISM_CLASS': open(join(uniprot_meta_dp, 'uniprot_organism_class.txt'), 'w'),
        'OTHER_ID': open(join(uniprot_meta_dp, 'uniprot_other_ids.txt'), 'w'),
        'RELATED_KEYWORD': open(join(uniprot_meta_dp, 'uniprot_related_keywords.txt'), 'w'),
        'RELATED_PUBMED_ID': open(join(uniprot_meta_dp, 'uniprot_related_pubmed_ids.txt'), 'w'),
        'SPECIES': open(join(uniprot_meta_dp, 'uniprot_species.txt'), 'w'),
    }
    with open(join(data_root, 'uniprot', 'uniprot_metadata.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')

            # Fail if metadata type is not in the map
            if p not in meta_output_files:
                raise Exception(f'Predicate not recognized {p}')

            meta_output_files[p].write(line)
    for fd in meta_output_files.values():
        fd.close()


def write_drugbank_metadata():
    drugbank_meta_dp = join(meta_root, 'drug')
    makedirs(drugbank_meta_dp) if not isdir(drugbank_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(drugbank_meta_dp, 'drugbank_name.txt'), 'w'),
        'SYNONYM': open(join(drugbank_meta_dp, 'drugbank_synonym.txt'), 'w'),
        'TYPE': open(join(drugbank_meta_dp, 'drugbank_type.txt'), 'w'),
        'PRODUCT': open(join(drugbank_meta_dp, 'drugbank_product.txt'), 'w'),
        'PUBMED_ARTICLE': open(join(drugbank_meta_dp, 'drugbank_related_pubmed_ids.txt'), 'w'),
        'DIRECT_PARENT': open(join(drugbank_meta_dp, 'drugbank_direct_parent.txt'), 'w'),
        'KINGDOM': open(join(drugbank_meta_dp, 'drugbank_kingdom.txt'), 'w'),
        'SUPERCLASS': open(join(drugbank_meta_dp, 'drugbank_superclass.txt'), 'w'),
        'CLASS': open(join(drugbank_meta_dp, 'drugbank_class.txt'), 'w'),
        'SUBCLASS': open(join(drugbank_meta_dp, 'drugbank_subclass.txt'), 'w'),
        'ALTERNATIVE_PARENT': open(join(drugbank_meta_dp, 'drugbank_alternative_parent.txt'), 'w'),
        'SUBSTITUENT': open(join(drugbank_meta_dp, 'drugbank_substituent.txt'), 'w'),
        'PRODUCT_STAGE': open(join(drugbank_meta_dp, 'drugbank_product_stage.txt'), 'w')
    }
    files = [
        join(data_root, 'drugbank', 'db_meta.txt'),
        join(data_root, 'drugbank', 'db_product_stage.txt'),
        join(data_root, 'drugbank', 'db_classification.txt'),
    ]
    for f in files:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o = line.strip().split('\t')

                # Fail if metadata type is not in the map
                if p not in meta_output_files:
                    raise Exception(f'Predicate not recognized {p}')

                meta_output_files[p].write(line)

    for fd in meta_output_files.values():
        fd.close()

    with open(join(data_root, 'drugbank', 'db_pathways.txt'), 'r') as fd:
        with open(join(pathway_properties_root, 'pathway_category.txt'), 'w') as output:
            for line in fd:
                s, p, o = line.strip().split('\t')
                if p == 'PATHWAY_CATEGORY':
                    output.write(line)


def write_pathway_go_annotations():
    pred_file_map = {
        'GO_BP': open(join(pathway_properties_root, 'pathway_go_biological_processes.txt'), 'w'),
        'GO_CC': open(join(pathway_properties_root, 'pathway_go_cellular_components.txt'), 'w'),
        'GO_MF': open(join(pathway_properties_root, 'pathway_go_molecular_functions.txt'), 'w')
    }

    with open(join(data_root, 'reactome', 'reactome_go_mapping.txt'), 'r') as fd:
        for line in fd:
            protein, map_type, goid, pathway_id, species = line.strip().split('\t')
            pred_file_map[map_type].write(f'{pathway_id}\tPATHWAY_{map_type}\t{goid}\n')
    
    for f in pred_file_map.values():
        f.close()


def write_mesh_metadata():
    mesh_meta_dp = join(meta_root, 'disease')
    makedirs(mesh_meta_dp) if not isdir(mesh_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(mesh_meta_dp, 'mesh_name.txt'), 'w'),
        'TYPE': open(join(mesh_meta_dp, 'mesh_type.txt'), 'w')
    }
    files = [
        join(data_root, 'mesh', 'mesh_disease_meta.txt'),
        join(data_root, 'mesh', 'mesh_scr_disease_meta.txt')
    ]
    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                s, p, o = line.strip().split('\t')

                # Fail if metadata type is not in the map
                if p not in meta_output_files:
                    raise Exception(f'Predicate not recognized {p}')

                meta_output_files[p].write(line)

    for fd in meta_output_files.values():
        fd.close()


def copy_folder(src_folder, dst_folder, included_files=[]):
    for f in listdir(src_folder):
        if f in included_files:
            copy(join(src_folder, f), dst_folder)


def compress_folder(folder):
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.gz'):
                continue
            src_fp = join(root, fp)
            dst_fp = join(root, fp+'.gz')
            with open(src_fp, 'rb') as f_in, gzip.open(dst_fp, 'wb') as f_out:
                f_out.writelines(f_in)
            remove(src_fp)


def generate_core_links():
    output = open(join(core_root, 'biokg.links.tsv'), 'w')
    for f in listdir(links_root):
        if f.endswith('.txt') and f != 'README.txt':
            with open(join(links_root, f), 'r') as fd:
                for line in fd:
                    output.write(line)
    output.close()

def generate_props(name, folder):
    output = open(join(core_root, f'biokg.properties.{name}.tsv'), 'w')
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.txt') and fp != 'README.txt':
                with open(join(root, fp)) as fd:
                    for line in fd:
                        output.write(line)
    output.close()

def generate_core_props():
    generate_props('protein', protein_properties_root)
    generate_props('drug', drug_properties_root)
    generate_props('cell', cell_properties_root)
    generate_props('pathway', pathway_properties_root)
    generate_props('disease', disease_properties_root)
    generate_props('genetic_disorder', join(properties_root, 'genetic_disorders'))



def generate_meta(name, folder):
    output = open(join(core_root, f'biokg.metadata.{name}.tsv'), 'w')
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.txt') and fp != 'README.txt':
                with open(join(root, fp)) as fd:
                    for line in fd:
                        output.write(line)
    output.close()


def generate_core_metadata():
    generate_meta('protein', join(meta_root, 'protein'))
    generate_meta('drug', join(meta_root, 'drug'))
    generate_meta('disease', join(meta_root, 'disease'))
    generate_meta('pathway', join(meta_root, 'pathway'))


def compile_graph():
    # Get the set of proteins, drugs and diseases to use and write metadata
    print('Writing Metadata')
    protein_set = get_all_proteins()
    write_uniprot_metadata()

    drug_set = get_all_drugs()
    write_drugbank_metadata()

    disease_set = get_all_mesh_diseases()
    write_mesh_metadata()

    genetic_disorders = get_all_genetic_disorders()
    # Get Links
    print('Writing Links')
    triples = get_all_unique_ppi(protein_set)
    print(f'{len(triples)} protein protein interactions')
    write_triples(triples, join(links_root, 'ppi.txt'))

    write_ppi_by_species()
    get_all_protein_sequence_annotations(protein_set)

    triples, function_triples = get_all_protein_drug_interactions(protein_set, drug_set)
    print(f'{len(triples)} drug targets')
    write_triples(triples, join(links_root, 'dpi.txt'))
    write_triples(function_triples, join(links_root, 'drug_protein_function.txt'))

    triples = get_all_protein_pathway_associations(protein_set)
    print(f'{len(triples)} protein pathway associations')
    write_triples(triples, join(links_root, 'protein_pathway.txt'))

    triples = get_all_protein_disease_associations(protein_set, disease_set)
    print(f'{len(triples)} protein disease associations')
    write_triples(triples, join(links_root, 'protein_disease.txt'))

    triples, tissue_structure = get_all_protein_expressions(protein_set)
    print(f'{len(triples)} protein tissue associations')
    write_triples(triples, join(protein_properties_root, 'protein_expression.txt'))
    write_triples(tissue_structure, join(cell_properties_root, 'cell_tissue_membership.txt'))

    triples = get_protein_complex_rels(protein_set)
    print(f'{len(triples)} protein complex rels')
    write_triples(triples, join(links_root, 'protein_complex.txt'))

    triples, tl_triples = get_complex_pathway_rels()
    print(f'{len(triples)} complex pathway rels')
    write_triples(triples, join(links_root, 'complex_pathway.txt'))
    write_triples(tl_triples, join(links_root, 'complex_top_level_pathway.txt'))

    triples = get_pathway_rels()
    print(f'{len(triples)} pathway rels')
    reactome_meta = join(meta_root, 'pathway')
    makedirs(reactome_meta) if not isdir(reactome_meta) else None
    write_triples(triples, join(reactome_meta, 'pathway_parent.txt'))

    triples = get_all_drug_pathway_associations(drug_set)
    print(f'{len(triples)} drug pathway associations')
    write_triples(triples, join(links_root, 'drug_pathway.txt'))

    triples = get_all_drug_drug_interactions(drug_set)
    print(f'{len(triples)} drug drug associations')
    write_triples(triples, join(links_root, 'ddi.txt'))

    triples = get_all_drug_disease_associations(drug_set, disease_set)
    print(f'{len(triples)} drug disease associations')
    write_triples(triples, join(links_root, 'drug_disease.txt'))

    triples = get_all_disease_pathway_associations(disease_set)
    print(f'{len(triples)} disease pathway associations')
    write_triples(triples, join(links_root, 'disease_pathway.txt'))

    triples = get_all_disease_genetic_disorder(disease_set, genetic_disorders)
    print(f'{len(triples)} disease genetic disorders')
    write_triples(triples, join(links_root, 'disease_genetic_disorders.txt'))
    # Get Properties
    print('Writing properties')
    triples, indication_triples = get_all_drug_side_effects(drug_set)
    print(f'{len(triples)} drug side effect associations')
    write_triples(triples, join(drug_properties_root, 'drug_sideeffect.txt'))
    write_triples(indication_triples, join(drug_properties_root, 'drug_indications.txt'))

    triples = get_all_drug_atc_codes(drug_set)
    print(f'{len(triples)} drug atc codes')
    write_triples(triples, join(drug_properties_root, 'drug_atc_codes.txt'))

    triples = get_disease_tree(disease_set)
    print(f'{len(triples)} disease tree')
    write_triples(triples, join(disease_properties_root, 'disease_tree.txt'))

    write_protein_cellline_expressions(protein_set)
    write_pathway_go_annotations()

    copy(
        join(data_root, 'medgen', 'mim_categories.txt'),
        mim_properties_root
    )

    # Copy Other datasets
    print('Copying other files')
    other_ctd = join(other_root, 'ctd')
    makedirs(other_ctd) if not isdir(other_ctd) else None
    copy_folder(
        join(data_root, 'ctd'),
        other_ctd,
        included_files=[
            'ctd_disease_biological_process.txt',
            'ctd_disease_cellular_component.txt',
            'ctd_disease_molecular_function.txt',
            'ctd_drug_phenotype.txt',
            'ctd_drug_protein_interactions.txt',
            'ctd_disease_kegg_pathway_association.txt',
            'ctd_disease_reactome_pathway_association.txt'
        ]
    )

    other_cello = join(other_root, 'cellosaurus')
    makedirs(other_cello) if not isdir(other_cello) else None
    copy_folder(
        join(data_root, 'cellosaurus'),
        other_cello,
        included_files=[
            'cl_cat.txt',
            'cl_geo.txt',
            'cl_map.txt',
            'cl_pmid.txt'
        ]
    )

    other_drugbank = join(other_root, 'drugbank')
    makedirs(other_drugbank) if not isdir(other_drugbank) else None
    copy_folder(
        join(data_root, 'drugbank'),
        other_drugbank,
        included_files=[
            'db_mechanism_or_action.txt',
            'db_mesh.txt'
        ]
    )

    other_hpa = join(other_root, 'hpa')
    makedirs(other_hpa) if not isdir(other_hpa) else None
    copy_folder(
        join(data_root, 'hpa'),
        other_hpa,
        included_files=[
            'hpa_antibodies.txt'
        ]
    )

    other_kegg = join(other_root, 'kegg')
    makedirs(other_kegg) if not isdir(other_kegg) else None
    copy_folder(
        join(data_root, 'kegg'),
        other_kegg,
        included_files=[
            'gene_network.txt',
            'glycan_pathway.txt',
            'network_disease.txt',
            'network_drug.txt',
            'network_pathway.txt',
            'disease_meta.txt'
        ]
    )

    other_mesh = join(other_root, 'mesh')
    makedirs(other_mesh) if not isdir(other_mesh) else None
    copy_folder(
        join(data_root, 'mesh'),
        other_mesh,
        included_files=[
            'mesh_disease_concept_heading.txt',
            'mesh_drug_concept_heading.txt',
            'mesh_scr_drug_meta.txt',
            'mesh_drug_tree.txt',
            'mesh_drug_meta.txt'
        ]
    )

    other_phos = join(other_root, 'phosphosite')
    makedirs(other_phos) if not isdir(other_phos) else None
    copy_folder(
        join(data_root, 'phosphosite'),
        other_phos,
        included_files=[
            'phosphorylation_site.txt'
        ]
    )

    other_reactome = join(other_root, 'reactome')
    makedirs(other_reactome) if not isdir(other_reactome) else None
    copy_folder(
        join(data_root, 'reactome'),
        other_reactome,
        included_files=[
            'reactome_go_mapping.txt',
            'reactome_isoform_pathway.txt'
        ]
    )

    other_sider = join(other_root, 'sider')
    makedirs(other_sider) if not isdir(other_sider) else None
    copy_folder(
        join(data_root, 'sider'),
        other_sider,
        included_files=[
            'sider_indications_meta.txt',
            'sider_effects_meta.txt'
        ]
    )

    print('Extracting core graph')
    generate_core_links()
    generate_core_props()
    generate_core_metadata()

    # Gzip output
    print('Compressing output')
    compress_folder(output_root)

    #Copy readmes to output
    copy('links_description.txt', join(links_root, 'README.txt'))
    copy('properties_description.txt', join(properties_root, 'README.txt'))
    copy('meta_description.txt', join(meta_root, 'README.txt'))


if __name__ == '__main__':
    compile_graph()
