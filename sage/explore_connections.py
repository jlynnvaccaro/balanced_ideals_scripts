import utils
import sys
from collections import defaultdict


def len_distr(gens, lens):
    d = defaultdict(int)
    for g in gens:
        d[len(g)] += 1
    return tuple([d[l] for l in lens])

d = utils.load_json_data(sys.argv[1])
middle_row_len = d["max_len"]//2



hist = defaultdict(int)
hist_56 = defaultdict(int)
for b in d["balanced_ideals"]:
    max_len = 0
    min_len = d["max_len"]
    for g in b["gen"]:
        len_g = len(g)
        if len_g>max_len:
            max_len = len_g
        if len_g<min_len:
            min_len = len_g
    if min_len >= 5 and max_len==6:
        tuple_56 = len_distr(b["gen"],(5,6))
        hist_56[tuple_56] += 1
    hist[(min_len,max_len)]+=1
print(hist)
print(hist_56)
