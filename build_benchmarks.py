from os import makedirs
from os.path import join, isdir
import gzip
from biodblinker import GeneNameLinker

def export_triplets(triplets, filepath):
    """ Export triplets to file

    Parameters
    ----------
    triplets : iter
        Triplets
    filepath : str
        File path
    """
    fd = open(filepath, "w")
    for s, p, o in triplets:
        fd.write(f"{s}\t{p}\t{o}\n")
    fd.close()


def build_benchmarks(preprocessed_dp, output_dp):
    """ Build BioKG benchmarks

    Parameters
    ----------
    preprocessed_dp : str
        Preprocessed directory path
    output_dp : str
        Output directory path
    """

    benchmarks_dp = join(output_dp, "benchmarks")
    makedirs(benchmarks_dp) if not isdir(benchmarks_dp) else None

    # ================================================================================================
    # Build DDI benchmarks
    # ------------------------------------------------------------------------------------------------
    # Benchmark #1: Drug-drug interactions and their effects on body minerals
    # Benchmark #2: Drug-drug interactions and their effects on therapeutic efficacy
    # ------------------------------------------------------------------------------------------------
    ddi_desc_fp = join(preprocessed_dp, "drugbank", "db_ddi.txt")

    mineral_triplets = set()
    efficacy_triplets = set()

    mineral_effects = {"calcemia", "glycemia", "kalemia", "atremia"}
    efficacy_effects = {"efficacy"}

    for line in open(ddi_desc_fp):
        d1, _, d2, ddi_effect = line.strip().split("\t")

        is_mineral_effect = bool(sum([ef in ddi_effect for ef in mineral_effects]))
        is_efficacy_effect = bool(sum([ef in ddi_effect for ef in efficacy_effects]))

        if is_mineral_effect:
            mineral_triplets.add((d1, ddi_effect, d2))
        if is_efficacy_effect:
            efficacy_triplets.add((d1, ddi_effect, d2))

    ddi_chem_fp = join(benchmarks_dp, "ddi_minerals.tsv")
    ddi_efficacy_fp = join(benchmarks_dp, "ddi_efficacy.tsv")
    export_triplets(mineral_triplets, ddi_chem_fp)
    export_triplets(efficacy_triplets, ddi_efficacy_fp)

    # ================================================================================================
    # Build DPI benchmarks
    # ------------------------------------------------------------------------------------------------
    # Benchmark #1: Drug-protein interactions of FDA approved drugs
    # Benchmark #2: Drug effects of protein expression for FDA approved drugs
    # ------------------------------------------------------------------------------------------------
    drug_stage_fp = join(preprocessed_dp, "drugbank", "db_product_stage.txt")
    drug_targets_fp = join(output_dp, "links", "dpi.txt.gz")
    drug_targets_effects_fp = join(preprocessed_dp, "ctd", "ctd_drug_protein_interactions.txt")

    fda_approved_drugs = set()
    fda_dpi = set()
    for line in open(drug_stage_fp):
        dr, _, stage = line.strip().split("\t")
        if stage == "approved":
            fda_approved_drugs.add(dr)

    for line in gzip.open(drug_targets_fp, 'rt'):
        dr, _, pr = line.strip().split("\t")
        if dr in fda_approved_drugs:
            fda_dpi.add((dr, "DPI", pr))

    dpi_fda_fp = join(benchmarks_dp, "dpi_fda.tsv")
    export_triplets(fda_dpi, dpi_fda_fp)

    dpi_inc_exp = set()
    dpi_dec_exp = set()
    for line in open(drug_targets_effects_fp):
        dr, effect_txt, pr, _, _ = line.strip().split("\t")
        if effect_txt == "INCREASES_EXPRESSION":
            dpi_inc_exp.add((dr, pr))
        if effect_txt == "DECREASES_EXPRESSION":
            dpi_dec_exp.add((dr, pr))
    dpi_exp_fp = join(benchmarks_dp, "dep_fda_exp.tsv")
    conflict_dpi = set.intersection(dpi_inc_exp, dpi_dec_exp)
    dpi_inc_exp_triplets = [[dr, "inc_expr", pr] for dr, pr in dpi_inc_exp if (dr, pr) not in conflict_dpi]
    dpi_dec_exp_triplets = [[dr, "dec_expr", pr] for dr, pr in dpi_dec_exp if (dr, pr) not in conflict_dpi]
    dpi_exp_triplets = dpi_inc_exp_triplets + dpi_dec_exp_triplets
    dep_exp_fda_triplets = [[dr, eff, pr] for dr, eff, pr in dpi_exp_triplets if dr in fda_approved_drugs]
    export_triplets(dep_exp_fda_triplets, dpi_exp_fp)
    # ================================================================================================
    # Build Phosphorylation benchmark
    # ------------------------------------------------------------------------------------------------
    # Benchmark #1: Human Kinase Substrate phosphorylation interactions from phosphosite plus and Cutillas 20
    # ------------------------------------------------------------------------------------------------
    linker = GeneNameLinker()
    human_uniprot_ids = set()
    psp_set = set()
    cutillas_set = set()
    psp_phos_fp = join(preprocessed_dp, 'phosphosite', 'kinase_substrate.txt')
    cutillas_fp = join(preprocessed_dp, 'cutillas', 'phosphorylation.txt')
    phos_fp = join(benchmarks_dp, 'phosphorylation.tsv')
    with open(join(preprocessed_dp, 'uniprot', 'uniprot_metadata.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'SPECIES' and o == 'HUMAN':
                human_uniprot_ids.add(s)
    
    with open(psp_phos_fp, 'r') as fd:
        for line in fd:
            kin, _, sub = line.strip().split('\t')[:3]
            if kin in human_uniprot_ids and sub in human_uniprot_ids:
                psp_set.add((kin, sub))

    with open(cutillas_fp, 'r') as fd:
        for line in fd:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                print(line)
                continue
            kin_gene, _, sub_gene = parts[:3]
            kin_uniprot_ids, sub_uniprot_ids = linker.convert_gene_name_to_uniprot([kin_gene, sub_gene])
            for kin in kin_uniprot_ids:
                if kin in human_uniprot_ids:
                    for sub in sub_uniprot_ids:
                        if sub in human_uniprot_ids:
                            cutillas_set.add((kin, sub))

    phosphorylation_set = psp_set.union(cutillas_set)
    phos_triplets = [[kin, 'phosphorylates', sub] for kin, sub in phosphorylation_set]
    export_triplets(phos_triplets, phos_fp)


def main():
    preprocessed_dp = "./data/preprocessed"
    output_dp = "./data/output"
    build_benchmarks(preprocessed_dp, output_dp)


if __name__ == '__main__':
    main()
