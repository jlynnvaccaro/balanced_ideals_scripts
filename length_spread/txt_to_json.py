import utils
import sys
import os

def write_ideals_to_json(infile, outfile):
    basename = os.path.basename(infile)
    name = os.path.splitext(basename)
    d = {"name":name,"summands":[],"ideals":[]}
    with open(infile, "rt") as f:
        summands = f.readline().strip().split("x")
        d["summands"] = [s.strip() for s in summands]
        for line in f:
            line=line.strip()
            if "gen:" in line:
                d["ideals"].append([g for g in utils.get_gen_set(line)])
    utils.save_json_data(outfile, d)
    print("Success! Wrote {} ideals to {}".format(infile,outfile))

if __name__=="__main__":
    if len(sys.argv) < 3:
        raise TypeError("Requires 3 input arguments")
    infile = sys.argv[1]
    outfile = sys.argv[2]
    write_ideals_to_json(infile, outfile)




