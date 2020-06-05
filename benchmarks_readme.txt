ddi_minerals.tsv
This benchmark contains drug drug interactions and the effect they have on minerals
Format
<drug> <mineral_effect> <drug>
Undirected

ddi_efficacy.tsv
This benchmark contains drug drug interactions and the effect they have on theraputic efficacy
Format
<drug> <efficact_effect> <drug>
Undirected

dpi_fda.tsv
This benchmark contains drug protein interactions for FDA approved drugs
Format
<drug> DPI <protein>
Directed

dep_fda_exp.tsv
This benchmark contains drug protein interaction and the effect they have on protein expression for FDA approved drugs
Format
<drug> inc_expr|dec_expr <protein>
Directed

phosphorylation.tsv
This benchmark contains protein protein phosphorylation interactions
Format
<kinase_protein> phosphorylates <substrate_protein> <substrate_site>
Directed