from biolink import KEGGLinker, SiderLinker
from os import makedirs
from os.path import join, isdir

kegg_linker = KEGGLinker()
sider_linker = SiderLinker()
data_root = 'data/preprocessed'
output_root = 'data/biokg'
makedirs(output_root) if not isdir(output_root) else None


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
    mesh_meta = join(data_root, 'mesh', 'mesh_metadata.txt')

    mesh_diseases = set()

    with open(mesh_meta, 'r') as fd:
        for line in fd:
            disease, meta, value = line.strip().split('\t')
            if meta == 'TYPE' and value == 'DISEASE':
                mesh_diseases.add(disease)

    return mesh_diseases

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
            unique_triples.append((s, 'PRO_PRO_INT', o))

    return unique_triples


def get_all_unique_phosphorylations(protein_set):
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
    ks_fp = join(data_root, 'phosphosite', 'kinase_substrate.txt')
    ppis = set()

    with open(ks_fp, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            ppis.add((s, o))

    unique_triples = []
    for s, o in ppis:
        if s in protein_set and o in protein_set:
            unique_triples.append((s, 'PRO_PRO_PHOS', o))

    return unique_triples


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
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    pdis = set()
    carriers = set()
    enzymes = set()
    transporters = set()
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

    with open(ctd_drug_protein, 'r') as fd:
        for line in fd:
            s, p, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                pdis.add((o, s))

    unique_triples = []
    for protein, drug in pdis:
        if protein in protein_set and drug in drug_set:
            unique_triples.append((protein, 'PRO_DRG_TGT', drug))

    unique_carriers = []
    for protein, drug in carriers:
        if protein in protein_set and drug in drug_set:
            unique_carriers.append((protein, 'PRT_DRG_CAR', drug))
    
    unique_enzymes = []
    for protein, drug in enzymes:
        if protein in protein_set and drug in drug_set:
            unique_enzymes.append((protein, 'PRT_DRG_ENZ', drug))

    unique_transporters = []
    for protein, drug in transporters:
        if protein in protein_set and drug in drug_set:
            unique_transporters.append((protein, 'PRT_DRG_TRANS', drug))
    
    return unique_triples, unique_carriers, unique_enzymes, unique_transporters


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

    # with open(uniprot_facts, 'r') as fd:
    #     for line in fd:
    #         s, p, o = line.strip().split('\t')
    #         if p == 'RELATED_DISEASE':
    #             o = o.split(':')[1]
    #             pdis.add((s, o))

    with open(ctd_protein_disease, 'r') as fd:
        for line in fd:
            s, p, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                pdis.add((s, o))

    unique_triples = []
    for protein, disease in pdis:
        if protein in protein_set and disease in disease_set:
            unique_triples.append((protein, 'PRO_DIS_ASS', disease))

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
    reactome_iso_pathway = join(data_root, 'reactome', 'reactome_isoform_pathway.txt')
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

    with open(reactome_iso_pathway, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')

            ppis.add((s.split('-')[0], o))

    unique_triples = []
    for protein, pathway in ppis:
        if protein in protein_set:
            unique_triples.append((protein, 'PRO_PTH_ASS', pathway))

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
    db_interactions = join(data_root, 'drugbank', 'db_interactions.txt')
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
            unique_triples.append((d1, 'DRG_DRG_INT', d2))

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
        join(data_root, 'ctd', 'ctd_drug_reactome_pathway_association.txt')
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
            unique_triples.append((drug, 'DRG_PTH_ASS', pathway))

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
            unique_triples.append((drug, 'DRG_DIS_ASS', disease))
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
            unique_triples.append((drug, 'DRG_SE_ASS', effect))

    for drug, indication in indications:
        if drug in drug_set:
            unique_triples.append((drug, 'DRG_IND_ASS', indication))
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
            unique_triples.append((disease, 'DIS_PTH_ASS', pathway))
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
    with open(hpa_tissue_expression, 'r') as fd:
        for line in fd:
            p, t, level = line.strip().split('\t')
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
    for protein, tissue in low_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'EXP_LOW', tissue))
    for protein, tissue in medium_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'EXP_MED', tissue))
    for protein, tissue in high_expression:
        if protein in protein_set:
            unique_triples.append((protein, 'EXP_HIGH', tissue))
    return unique_triples


def write_triples(triples, output_fp):
    with open(output_fp, 'w') as output:
        for s, p, o in triples:
            output.write(f'{s}\t{p}\t{o}\n')

def filter_ctd_drug_protein(protein_set):
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    ctd_drug_protein_filtered = join(output_root, 'ctd_drug_protein_interactions_SWISSPORT.txt' )
    with open(ctd_drug_protein_filtered, 'w') as output_fd:
        with open(ctd_drug_protein, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')[:4]
                if o in protein_set:
                    output_fd.write(line)

protein_set = get_all_proteins()
#protein_set = get_proteins_by_metadata('SPECIES', ['HUMAN'])

drug_set = get_all_drugs()
disease_set = get_all_mesh_diseases()

triples = get_all_unique_ppi(protein_set)
print(f'{len(triples)} protein protein interactions')
write_triples(triples, join(output_root, 'protein_protein.txt'))

triples = get_all_unique_phosphorylations(protein_set)
print(f'{len(triples)} protein protein phosphorylations')
write_triples(triples, join(output_root, 'protein_protein_phos.txt'))

triples, cars, enzs, trans = get_all_protein_drug_interactions(protein_set, drug_set)
print(f'{len(triples)} protein drug targets')
write_triples(triples, join(output_root, 'protein_drug.txt'))
write_triples(cars, join(output_root, 'protein_drug_carriers.txt'))
write_triples(enzs, join(output_root, 'protein_drug_enzymes.txt'))
write_triples(trans, join(output_root, 'protein_drug_transporters.txt'))

triples = get_all_protein_pathway_associations(protein_set)
print(f'{len(triples)} protein pathway associations')
write_triples(triples, join(output_root, 'protein_pathway.txt'))

triples = get_all_protein_disease_associations(protein_set, disease_set)
print(f'{len(triples)} protein disease associations')
write_triples(triples, join(output_root, 'protein_disease.txt'))

triples = get_all_protein_expressions(protein_set)
print(f'{len(triples)} protein tissue associations')
write_triples(triples, join(output_root, 'protein_expression.txt'))

triples = get_all_drug_pathway_associations(drug_set)
print(f'{len(triples)} drug pathway associations')
write_triples(triples, join(output_root, 'drug_pathway.txt'))

triples = get_all_drug_drug_interactions(drug_set)
print(f'{len(triples)} drug drug associations')
write_triples(triples, join(output_root, 'drug_drug.txt'))

triples = get_all_drug_disease_associations(drug_set, disease_set)
print(f'{len(triples)} drug disease associations')
write_triples(triples, join(output_root, 'drug_disease.txt'))

triples = get_all_drug_side_effects(drug_set)
print(f'{len(triples)} drug side effect associations')
write_triples(triples, join(output_root, 'drug_sideeffect.txt'))

triples = get_all_disease_pathway_associations(disease_set)
print(f'{len(triples)} disease pathway associations')
write_triples(triples, join(output_root, 'disease_pathway.txt'))
