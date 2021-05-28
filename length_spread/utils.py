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