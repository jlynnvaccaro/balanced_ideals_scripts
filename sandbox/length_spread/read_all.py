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
    if "summary" in datafile or "tail" in datafile or "png" in datafile:
        print("Skipping",datafile)
        continue
    # if datafile != "A1A2.txt":
        # continue
    print(datafile)

    if datafile[-4:] ==".txt":
        f = open(datafile, "rt")
    else:
        f = gzip.open(datafile, 'rt')

    for line in f:
        line=line.strip()
        if "x" in line:
            print("----",line)

