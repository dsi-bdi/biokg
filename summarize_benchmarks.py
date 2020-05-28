from os.path import join
from os.path import basename
import numpy as np


def summarize_benchmark(filepath):
    fd = open(filepath)
    benchmarks_data = [line.strip().split("\t") for line in fd]
    subs = set()
    rels = set()
    objs = set()
    pair = set()

    for s, p, o in benchmarks_data:
        subs.add(s)
        rels.add(p)
        objs.add(o)
        pair.add((s, o)) if s > o else pair.add((o, s))
    ents = set.union(subs, objs)
    nb_subs = len(subs)
    nb_objs = len(objs)
    nb_rels = len(rels)
    nb_ents = len(ents)
    nb_pairs = len(pair)
    nb_tps = len(benchmarks_data)
    return {"nb_ents": nb_ents, "nb_rels": nb_rels, "nb_subs": nb_subs, "nb_objs": nb_objs, "nb_triplets": nb_tps, "nb_pairs": nb_pairs}


if __name__ == '__main__':
    data_dp = "./data"
    benchmarks_dp = join(data_dp, "output", "benchmarks")

    ddi_efficacy_fp = join(benchmarks_dp, "ddi_efficacy.tsv")
    ddi_minerals_fp = join(benchmarks_dp, "ddi_minerals.tsv")
    dep_fda_exp_fp = join(benchmarks_dp, "dep_fda_exp.tsv")
    dpi_fda_fp = join(benchmarks_dp, "dpi_fda.tsv")

    benchmark_filepath_list = [ddi_efficacy_fp, ddi_minerals_fp, dep_fda_exp_fp, dpi_fda_fp]
    for benchmark_filepath in benchmark_filepath_list:
        summary = summarize_benchmark(benchmark_filepath)
        print("======================================================================")
        print(f"= {basename(benchmark_filepath)}")
        print("=================================")
        print(f"= # Triplets:  {summary['nb_triplets']}")
        print(f"= # Entities:  {summary['nb_ents']}")
        print(f"= # Relations: {summary['nb_rels']}")
        print(f"= # Subjects:  {summary['nb_subs']}")
        print(f"= # Objects:   {summary['nb_objs']}")
        print(f"= # SO Pairs:   {summary['nb_pairs']}")
