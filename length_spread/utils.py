def get_ideal_stats(line):
    L = line[28:].strip().split()
    m = 0
    for gen in L:
        gen_l = len(gen)
        if gen_l>m:
            m = gen_l
    return m, len(L)