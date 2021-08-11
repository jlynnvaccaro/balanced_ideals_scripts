import json
import datetime

def add_header(d,version="0.0.1"):
    """Adds the json header to files"""
    now = datetime.datetime.now()
    json_dt = now.isoformat()
    d["timestamp"] = json_dt
    d["creator"] = __file__
    d["version"] = version

def get_ideal_stats(line):
    """Gets the max_length of the ideal and the number of generators"""
    L = line.strip().split(":")
    L = L[-1].strip().split()
    m = 0
    for gen in L:
        gen_l = len(gen)
        if gen_l>m:
            m = gen_l
    return m, len(L)

def get_gen_set(line):
    L = line.strip().split(":")
    L = L[-1].strip().split()
    return L

def load_json_data(filename):
    with open(filename, "rt") as f:
        return json.load(f)

def save_json_data(filename, data):
    with open(filename, "wt") as f:
        json.dump(data, f)

def simple_list():
        return ["A1.txt","A2.txt","A3.txt","A4.txt",
        "B2.txt","B3.txt","D4.txt","G2.txt"]

def fast_data_list():
    return ["A1A1A1A1A1A1A1.txt",
    "A1A1A2.txt",
    "A2B2.txt",
    "A1A1A1A1A1A1.txt",
    "A1A1A3.txt",
    "A1B2B2.txt",
    "B2B2.txt",
    "A1A1A1A1A1.txt",
    "A1A1B2.txt",
    "A1B2G2.txt",
    "A2G2.txt",
    "A1A1A1A1A2.txt",
    "A1A1G2.txt",
    "A1B2.txt",
    "B2G2.txt",
    "D4.txt",
    "A1A1A1A1B2.txt",
    "A1A1.txt",
    "A1B3.txt",
    "A2.txt",
    "G2G2.txt",
    "A1A1A1A1.txt",
    "A1A2A2.txt",
    "A3B2.txt",
    "B2.txt",
    "A1A1A1A2.txt",
    "G2.txt",
    "A1A1A1B2.txt",
    "A1A2B2.txt",
    "A1G2.txt",
    "A1A1A1G2.txt",
    "A1A2G2.txt",
    "A1.txt",
    "B3.txt",
    "A1A1A1.txt",
    "A1A2.txt",
    "A2A2.txt",
    "A3.txt",
    "A1A1A2A2.txt",
    "A1A3.txt",
    "A2A3.txt",
    "A4.txt"]

def slow_data_list():
    return ["A1A2A3.txt.gz","A1A4.txt.gz","A1D4ideals.txt.gz","A2B3.txt.gz","A3G2.txt.gz","B2B3ideals.txt.gz"]