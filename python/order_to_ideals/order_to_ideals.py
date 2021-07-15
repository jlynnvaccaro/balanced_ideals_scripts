## 5/25/2021 J Vaccaro
## Script for calculating a linear regression between the order of a Weyl group and its number of balanced ideals.

import os

orders = []
ideals = []
DATAPATH = "/nas/share/ideals2021/data/"
os.chdir(DATAPATH)
for datafile in os.listdir("."):
    if datafile[-4:] != ".txt":
        continue
    print(datafile)
    order = 0
    ideal = 0
    with open(datafile, "rt") as f:
        for line in f:
            line=line.strip()
            if "Order:" in line:
                words = line.split()
                maybe_order = int(words[3])
                if maybe_order>order:
                    order = maybe_order
            elif "Found" in line:
                words = line.split()
                maybe_ideal = int(words[1])
                if maybe_ideal>ideal:
                    ideal = maybe_ideal
    orders.append(order)
    ideals.append(ideal)

import numpy