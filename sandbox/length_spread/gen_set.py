# Is the number of generators of a lower-half balanced ideal the maximum # of generators among all balanced ideals?
import utils
from collections import defaultdict
import gzip
import os
# datafile = "/nas/share/ideals2021/data/A1A1A1A1A1.txt"

big_list = []

# DATAPATH = "/home/jvacca4/data/"
DATAPATH = "/nas/share/ideals2021/data/"
os.chdir(DATAPATH)
for datafile in os.listdir("."):
    if "summary" in datafile or "tail" in datafile:
        print("Skipping",datafile)
        continue
    # if datafile != "A1A2.txt":
        # continue
    print(datafile)
    d = {"name":datafile,"summand":1,"order":0,"l(w0)":0,"num_gens":defaultdict(int),"Rh num_gens":defaultdict(int)}
    gen_set = set()


    if datafile[-4:] ==".txt":
        f = open(datafile, "rt")
    else:
        f = gzip.open(datafile, 'rt')

    for line in f:
        line=line.strip()
        if "x" in line:
            d["summand"] += line.count("x")
        if "Order:" in line:
            words = line.split()
            maybe_order = int(words[3])
            if maybe_order>d["order"]:
                d["order"] = maybe_order
            d["l(w0)"] = int(words[6])
        # elif "Found" in line:
        #     words = line.split()
        #     maybe_total = int(words[1])
        #     if maybe_total>d["total"]:
        #         d["total"] = maybe_total
        elif "gen:" in line:
            # for i in utils.get_gen_set(line):
            #     gen_set.add(i)
            m,n = utils.get_ideal_stats(line)
            d["num_gens"][n] += 1
            if m<=d["l(w0)"]/2:
                d["Rh num_gens"][n] += 1
    L = list(gen_set)
    d["gen_set"] = L
    big_list.append(d)


for grp in big_list:
    print(grp["name"])
    # Max # of generators in an ideal, and how many
    k1 = max(grp["num_gens"].keys())
    # print(k1, grp["num_gens"][k1])
    # # of 
    if len(grp["Rh num_gens"].keys())>0:
        k2 = max(grp["Rh num_gens"].keys())
        # print(k2, grp["num_gens"][k2])
        if k1==k2 and grp["num_gens"][k1] == grp["num_gens"][k2]:
            print("TRUE!")
        else:
            print("-------False.")
    else:
        print("No lower half ideals")