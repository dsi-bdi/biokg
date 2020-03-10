import re
from collections import defaultdict

pat_1 = re.compile('The risk or severity of .* can be increased when .* is combined with .*')
pat_2 = re.compile('The risk or severity of .* can be decreased when .* is combined with .*')

pat_3 = re.compile('.* may decrease effectiveness of .* as a diagnostic agent.')

pat_4 = re.compile('The therapeutic efficacy of .* can be increased when used in combination with .*')
pat_5 = re.compile('The therapeutic efficacy of .* can be decreased when used in combination with .*')

pat_6 = re.compile('.* may increase the .* activities of .*')
pat_7 = re.compile('.* may decrease the .* activities of .*')

pat_8 = re.compile('The serum concentration of .* can be increased when it is combined with .*')
pat_9 = re.compile('The serum concentration of .* can be decreased when it is combined with .*')

pat_10 = re.compile('The metabolism of .* can be increased when combined with .*')
pat_11 = re.compile('The metabolism of .* can be decreased when combined with .*')

pat_12 = re.compile('The excretion of .* can be decreased when combined with .*')
pat_13 = re.compile('The excretion of .* can be increased when combined with .*')

pat_14 = re.compile('.* can cause a decrease in the absorption of .* resulting in a reduced serum concentration and potentially a decrease in efficacy.')

pat_15 = re.compile('.* may decrease the excretion rate of .* which could result in a higher serum level.')
pat_16 = re.compile('.* may increase the excretion rate of .* which could result in a lower serum level and potentially a reduction in efficacy.')

pat_17 = re.compile('The protein binding of .* can be decreased when combined with .*')

pat_18 = re.compile('.* can cause an increase in the absorption of .* resulting in an increased serum concentration and potentially a worsening of adverse effects.')

pat_19 = re.compile('The risk of a hypersensitivity reaction to .* is increased when it is combined with .*')

pat_20 = re.compile('The bioavailability of .* can be increased when combined with .*')
pat_21 = re.compile('The bioavailability of .* can be decreased when combined with .*')

pat_22 = re.compile('The absorption of .* can be decreased when combined with .*')

pat_23 = re.compile('The serum concentration of the active metabolites of .* can be increased when .* is used in combination with .*')
pat_24 = re.compile('The serum concentration of the active metabolites of .* can be reduced when .* is used in combination with .*')
pat_25 = re.compile('The serum concentration of the active metabolites of .* can be decreased when .* is used in combination with .*')

pat_26 = re.compile('The serum concentration of .* an active metabolite of .* can be increased when used in combination with .*')
pat_27 = re.compile('The serum concentration of .* an active metabolite of .* can be decreased when used in combination with .*')

patterns = [
    pat_1, pat_2, pat_3, pat_4, pat_5,
    pat_6, pat_7, pat_8, pat_9, pat_10,
    pat_11, pat_12, pat_13, pat_14, pat_15,
    pat_16, pat_17, pat_18, pat_19, pat_20,
    pat_21, pat_22, pat_23, pat_24, pat_25,
    pat_26, pat_27
]

pat_hits = defaultdict(int)
with open('data/preprocessed/db_interactions.txt', 'r') as input_fd:
    for line in input_fd:
        parts = line.strip().split('\t')
        desc = parts[-1]
        pat_found = False
        for index, pattern in enumerate(patterns):
            if re.fullmatch(pattern, desc) is not None:
                pat_hits[index] += 1
                pat_found = True
                break

        if not pat_found:
            print(line.strip())

total_hits = 0
for index in range(len(patterns)):
    print(f'{index}\t{pat_hits[index]}\t{patterns[index]}')
    total_hits += pat_hits[index]

print(total_hits)