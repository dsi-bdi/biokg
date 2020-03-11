import re
from collections import defaultdict
from tqdm import tqdm


def sanatize_se_txt(txt):
    return txt.strip().replace(" ", "_").lower()


pat_1 = re.compile('The risk or severity of (?P<se>.*) can be (?P<mode>\S+)d when .* is combined with .*')
pat_2 = re.compile('.* may (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) of .* as a diagnostic agent.')
pat_3 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when used in combination with .*')
pat_4 = re.compile('.* may (?P<mode>\S+) the (?P<se>.*) of .*')
pat_5 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when it is combined with .*')
pat_6 = re.compile('.* can cause a decrease in the absorption of .* resulting in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) and potentially a decrease in efficacy.')
pat_7 = re.compile('.* may decrease the excretion rate of .* which could result in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*).')
pat_8 = re.compile('.* may increase the excretion rate of .* which could result in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) and potentially a reduction in efficacy.')
pat_9 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when combined with .*')
pat_10 = re.compile('.* can cause an increase in the absorption of .* resulting in an (?P<mode>\S+)d (?P<se>\S+\s?\w*\s?\w*) and potentially a worsening of adverse effects.')
pat_11 = re.compile('The risk of a (?P<se>\S+\s?\w*\s?\w*) to .* is (?P<mode>\S+)d when it is combined with .*')
pat_12 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when combined with .*')
pat_13 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of the active metabolites of .* can be (?P<mode>\S+)d when .* is used in combination with .*')
pat_14 = re.compile('The (?P<se>\S+\s?\w*\s?\w*) of .*, an active metabolite of .* can be (?P<mode>\S+)d when used in combination with .*')
pat_15 = re.compile('.* may (?P<mode>\S+) the central nervous system depressant (?P<se>\S+\s?\S*\s?\S*) of .*')

patterns = [
    pat_1, pat_2, pat_3, pat_4, pat_5,
    pat_6, pat_7, pat_8, pat_9, pat_10,
    pat_11, pat_12, pat_13, pat_14, pat_15
]

mode_map = {
    'reduced': "decrease",
    'increase': "increase",
    'higher': "increase",
    'decrease': "decrease",
    'reduce': "decrease",
    'lower': "decrease"
}

side_effect_name_map = {
    "central_nervous_system_depressant_(cns_depressant)_activities": 'cns_depression_activities',
    "cns_depression": 'cns_depression_activities',
    "cardiotoxic_activities": 'cardiotoxicity',
    "constipating_activities": 'constipation',
    "excretion": 'excretion_rate',
    "hyperkalemic_activities": 'hyperkalemia',
    "hypertensive_activities": 'hypertension',
}

pat_hits = defaultdict(int)
se_map = dict()
missed_lines = 0
missed_lines_list = []
mode_list = set()
out_fd = open('data/preprocessed/db_ddi_se.txt', 'w')
with open('data/preprocessed/db_interactions.txt', 'r') as input_fd:
    for line in tqdm(input_fd):
        d1, _, d2, desc = line.strip().split('\t')
        pat_found = False
        for pattern_index, pattern in enumerate(patterns):
            pg = re.match(pattern, desc)
            if pg is not None:
                se_name_list = []
                se_name = pg.group("se").lower()
                mode = pg.group("mode")
                has_word_activities = ("activities" in se_name)
                if has_word_activities:
                    se_name = se_name.replace(" activities", "")
                mode_name = mode_map[mode]
                if ", and" in se_name:
                    se_name_list = [sanatize_se_txt(se) for se in se_name.replace("and", "").split(", ")]
                elif "and" in se_name:
                    se_name_list = [sanatize_se_txt(se) for se in se_name.split(" and ")]
                else:
                    se_name_list = [sanatize_se_txt(se_name)]

                if has_word_activities:
                    se_name_list = [txt+"_activities" for txt in se_name_list]

                for side_effect in se_name_list:
                    if side_effect in side_effect_name_map:
                        side_effect = side_effect_name_map[side_effect]

                    out_fd.write(f"{d1}\t{d2}\t{mode_name}_{side_effect}\n")
                    if side_effect not in se_map:
                        se_map[side_effect] = 1
                    else:
                        se_map[side_effect] += 1

                pat_hits[pattern_index] += 1
                pat_found = True
                break

        if not pat_found:
            missed_lines += 1
            # print(line.strip())

out_fd.close()
total_hits = 0

for index in range(len(patterns)):
    # print(f'{index}\t{pat_hits[index]}\t{patterns[index]}')
    total_hits += pat_hits[index]

print(f'total hits: {total_hits} - missed lines:{missed_lines}')