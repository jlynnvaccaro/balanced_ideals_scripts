import utils
import sys
from collections import defaultdict

def bruhat_len(s):
    """Gives the length of the word"""
    if s=="1":
        return 0
    return len(s)

def print_hist(max_len, d):
    """Given a defaultdict, prints it in a nice order."""
    for i in range(max_len+1):
        print("{}: {}".format(i,d[i]))

def make_len_hist(L):
    """Make a histogram of lengths using defaultdict"""
    d = defaultdict(int)
    max_len = 0
    for x in L:
        x_len = bruhat_len(x)
        d[x_len] += 1
        if x_len > max_len:
            max_len = x_len
    return max_len, d

if __name__=="__main__":
    d = utils.load_json_data(sys.argv[1])
    max_len, hist = make_len_hist(d["core"])
    print("Type:",d["cartan_type"])
    print("Total size:",d["weyl_order"])
    print("Core size:",d["num_core"])
    print("Core length:",max_len)
    print("Core distribution:")
    print_hist(max_len,hist)
