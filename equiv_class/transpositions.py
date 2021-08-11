from collections import defaultdict
import equivalence_class
# from collections import defaultdict
class weyl_elt:
    def __init__(self,oneline):
        self.update_oneline(oneline)
    def update_oneline(self,oneline):
        self.oneline = oneline
        self.oneline_len = len(oneline)
        self.word = equivalence_class.to_word(oneline)
        if self.word=="1":
            self.word_len = 0
        else:
            self.word_len = len(self.word)
        self.find_triples()
    def __eq__(self,other):
        return self.word==other.word
    def __repr__(self):
        return "{}:{}".format(self.word,equivalence_class.matrix(self.oneline).sum().total_sum())
    def permute(self, other):
        new_oneline = [other.oneline[i] for i in self.oneline]
        return weyl_elt(new_oneline)


    def find_triples(self):
        num_triples = 0
        for i in range(self.oneline_len-2):
            if self.oneline[i]<self.oneline[i+1]<self.oneline[i+2]:
                num_triples+=1
        self.num_triples = num_triples

    def make_elts(self, d):
        for i in range(self.oneline_len-1):
            new_oneline = [x for x in self.oneline]
            new_oneline[i],new_oneline[i+1] = new_oneline[i+1],new_oneline[i]
            new_elt = weyl_elt(new_oneline)
            if new_elt.word_len < self.word_len:
                continue
            else:
                if new_elt not in d[new_elt.word_len]:
                    d[new_elt.word_len].append(new_elt)


def gen_all(n):
    D = {}
    root = weyl_elt(equivalence_class.to_oneline("1",n))
    D[0] = [root]
    i = 0
    while D[i]:
        D[i+1] = []
        print("----------------- start row",i)
        for elt in D[i]:
            print(elt)
            print(equivalence_class.matrix(elt.oneline).sum())
            elt.make_elts(D)
        print("----------------- end row",i)
        i += 1

from termcolor import colored
# import matplotlib.pyplot as plt
# import numpy as np
def print_matrices(elt, color="white"):
    """Prints out an element and some matrices"""
    w0_elt = elt.w0()
    s = "{} {}:{} {}".format(elt.word,elt.oneline,w0_elt.word,w0_elt.oneline)
    print(colored(s,color))
    print(colored("------",color))
    elt_sigma = elt.matrix_sigma()
    print(elt_sigma)
    w0_sigma = w0_elt.matrix_sigma()
    print(w0_sigma)
    print(colored("------",color))
    elt_delta = elt_sigma.dot()
    print(elt_delta)
    w0_delta = w0_sigma.dot()
    print(w0_delta)
    # plt.figure()
    # plt.imshow(np.array(elt_delta.data))
    # plt.imshow(np.array(w0_delta.data))
    # plt.show()
    # print(colored("-------------------------",color))
    # exit(0)


A = equivalence_class.equiv_tree(3)
A.make_equiv_class()
equiv = A.nodes("equiv")
print("\nEquivalence class")
print("Size:",len(equiv))
# for elt in equiv:
#     print_matrices(elt,"white")
# exit(0)

for elt in equiv:
    print(elt.word,end=", ")

A.make_ascending_equiv_class()
ascend = A.nodes("ascend")
print("\n\nEquivalence class (ascending only)")
print("Size:",len(ascend))
for elt in ascend:
    print(elt.word,end=", ")

A.make_symmetric_group()
sym = A.nodes("sym")
print("\n\nSymmetric group")
print("Size:",len(sym))
D = defaultdict(list)
for elt in sym:
    # print(elt.word,end=", ")
    D[elt.word_len].append(elt)
# print("")



for i in range(max(D.keys())//2+1):
    print(" Rank:",i)
    print("--------------------------------------------")
    for elt in D[i]:
        if elt in ascend:
            print_matrices(elt,"green")
        elif elt in equiv:
            print_matrices(elt,"blue")
        else:
            print_matrices(elt,"white")
    print("\n===========================================")

# hull_c = []
# # D = defaultdict(list)
# for elt in core:
#     # D[elt.word_len].append(elt)
#     w0_elt = elt.w0()
#     print(elt.word, elt.oneline, ":", w0_elt.word, w0_elt.oneline)
#     # print(elt.matrix_permutation())
#     elt_sigma = elt.matrix_sigma()
#     w0_sigma = w0_elt.matrix_sigma()
#     print(elt_sigma > w0_sigma)
#     print(elt_sigma < w0_sigma)
#     print(elt_sigma.total_sum(), w0_sigma.total_sum(), elt_sigma.total_sum()-w0_sigma.total_sum())
#     print(elt_sigma-w0_sigma)
#     print("---")
#     print(elt.matrix_sigma())
#     print("---")
#     print(w0_elt.matrix_sigma())
# # for i in range(max(D.keys())+1):
# #     for elt in D[i]:
# #         print(str(elt),end=" ")
# #     print("")

# # print("----------------")
# # gen_all(4)