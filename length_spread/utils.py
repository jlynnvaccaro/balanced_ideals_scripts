import json

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