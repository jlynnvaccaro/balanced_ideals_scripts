import utils
import sys

if len(sys.argv) < 3:
    raise TypeError("Requires 3 input arguments")

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "rt") as f:
    d = {"name":infile,"type":"","summands":[],"ideals":[]}
    with open(infile, "rt") as f:
        summands = f.readline().strip().split("x")
        d["type"] = summands[0]
        d["summands"] = [s for s in summands]
        for line in f:
            line=line.strip()
            if "gen:" in line:
                d["ideals"].append([g for g in utils.get_gen_set(line)])

utils.save_json_data(outfile, d)