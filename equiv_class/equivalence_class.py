import utils
global alphabet
alphabet = "abcdefghijklmnopqrstuvwxyz"

def to_oneline(s, n=4):
    """Converts a string from word rep to oneline rep in the Weyl group An"""
    L = [x for x in range(n+1)]
    if s == "1":
        return L
    for c in s[::-1]:
        i = alphabet.index(c)        
        L[i],L[i+1] = L[i+1],L[i]
    return L

def to_word(L):
    """Converts a oneline rep list to a word"""
    L_copy = [i for i in L]
    s=""
    for i in range(len(L)):
        ind = L_copy.index(i)
        while ind != i:
            s += alphabet[ind-1]
            L_copy[ind-1],L_copy[ind] = L_copy[ind],L_copy[ind-1]
            ind -= 1
    if s=="":
        return "1"
    return s

def is_increasing(L):
    """Checks whether a 3-elt numerical list is increasing"""
    if L[0]<L[1] and L[1]<L[2]:
        return True
    return False

def new_pattern_1(L,i=0):
    """Turns [1,2,3] into [2,1,3] starting at index i"""
    new_L = [x for x in L]
    new_L[i],new_L[i+1] = new_L[i+1],new_L[i]
    return new_L

def new_pattern_2(L,i=0):
    """Turns [1,2,3] into [1,3,2]"""
    return new_pattern_1(L,i+1)

class equiv_node:
    def __init__(self,oneline):
        self.oneline=oneline
        self.word = to_word(oneline)
        self.children = []
        self.confirmed_leaf = False
    def set_parent(self,parent):
        parent.add_child(self)
        self.parent=parent
    def add_child(self,child):
        self.children.append(child)
    def is_leaf(self):
        return not self.children
    def is_confirmed_leaf(self):
        return self.confirmed_leaf
    def get_children(self):
        return self.children
    def __eq__(self, other):
        return self.oneline==other.oneline
    def __repr__(self):
        return self.word

class equiv_tree:
    def __init__(self,n):
        self.root = equiv_node(to_oneline("1",n))
        self.range = n-1
    def equiv_nodes(self,L=None,elt=None):
        """Returns a list of equivalence class nodes"""
        if L is None:
            L = []
            return_L = True
        else: 
            return_L = False
        if elt is None:
            elt = self.root
        L.append(elt)
        for child in elt.get_children():
            self.equiv_nodes(L,child)
        if return_L:
            return L
    def equiv_leafs(self,L=None,elt=None):
        """Returns a list of equivalence class leafs"""
        if L is None:
            L=[]
            if elt is None:
                elt=self.root
            self.equiv_nodes(L,elt)
        return [x for x in L if x.is_leaf()]
    def make_equiv_class(self):
        """Generates the nodes related to the equivalence class"""
        num_loops = 0
        while True:
            num_loops += 1
            nodes = self.equiv_nodes()
            leafs = self.equiv_leafs(nodes)
            if all([x.is_confirmed_leaf() for x in leafs]):
                print("All done, found {} nodes in the equivalence class".format(len(nodes)))
                return num_loops
            for leaf in leafs:
                if leaf.is_confirmed_leaf():
                    continue
                for i in range(self.range):
                    if is_increasing(leaf.oneline[i:i+3]) or is_increasing(new_pattern_1(leaf.oneline[i:i+3])):
                        new_node = equiv_node(new_pattern_1(leaf.oneline,i))
                        if new_node not in nodes:
                            new_node.set_parent(leaf)
                            nodes.append(new_node)
                    if is_increasing(leaf.oneline[i:i+3]) or is_increasing(new_pattern_2(leaf.oneline[i:i+3])):
                        new_node = equiv_node(new_pattern_2(leaf.oneline,i))
                        if new_node not in nodes:
                            new_node.set_parent(leaf)
                            nodes.append(new_node)
                if leaf.is_leaf():
                    # print("Confirmed leaf:",leaf)
                    leaf.confirmed_leaf = True

def equiv_class(n=4, L=None, elt=None):
    """Recursively defines the equivalence class of any elt.
    Default is 1"""
    if elt is None:
        elt = to_oneline("1",n)
    if L is None:
        L = [elt]
    orig_len_L = len(L)
    for i in range(n-1):
        if is_increasing(elt[i:i+3]):
            new_elt_1 = new_pattern_1(elt,i)

            if new_elt_1 not in L:
                L.append(new_elt_1)
                equiv_class(n,L,new_elt_1)
            
            new_elt_2 = new_pattern_2(elt,i)
            if new_elt_2 not in L:
                L.append(new_elt_2)
                equiv_class(n,L,new_elt_2)
        if is_increasing(new_pattern_1(elt[i:i+3])):
            new_elt_1 = new_pattern_1(elt,i)

            if new_elt_1 not in L:
                L.append(new_elt_1)
                equiv_class(n,L,new_elt_1)
        if is_increasing(new_pattern_2(elt[i:i+3])):
            new_elt_2 = new_pattern_2(elt,i)

            if new_elt_2 not in L:
                L.append(new_elt_2)
                equiv_class(n,L,new_elt_2)
    if orig_len_L==1:
        return [to_word(x) for x in L]

# L_abc = to_oneline('abcdhdcadeded',10)
# print(L_abc)
# abc = to_word(L_abc)
# print(abc)

# for i in range(1,5):
#     L = equiv_class(i)
#     print(len(L))
A = equiv_tree(3)
print(A.make_equiv_class())
print(A.equiv_nodes())


A4_json_core_list= ["1","d","c","b","a","dc","cd","bd","cb","bc","ad","ac","ba","ab","cdc","cbd","bcb","dcb","bdc","adc","acd","bad","bac","abd","aba","acb","abc","dcbd","dbcb","acdc","badc","abad","acbd","abcb","abac"]
A4_json_core = {to_word(to_oneline(x,4)) for x in A4_json_core_list}
A4_equiv_core = set(["1", "a", "ba", "cba", "cbab", "bab", "babd", "cbac", "bac", "badc", "bad", "ac", "ad", "acd", "b", "cb", "acb", "ab", "abd", "acbd", "cbc", "bd", "c", "bc", "bdc", "dc", "adc", "adcd", "dcd", "d", "cd", "bcd", "cbcd", "cbd", "bdcd"])
A4_json_reverse = {to_word(to_oneline(x[::-1],4)) for x in A4_json_core}
print(A4_json_reverse==A4_equiv_core)

A5_json_core_list = ["1","e","d","c","b","a","ed","de","ce","dc","cd","be","bd","cb","bc","ae","ad","ac","ba","ab","ded","dce","cdc","edc","ced","bed","bde","cbe","cbd","bce","bcb","dcb","bdc","bcd","aed","ade","ace","adc","acd","bae","bad","bac","abe","abd","aba","acb","abc","edce","ecdc","dced","bded","cbed","bcbe","dcbe","dcbd","bdce","dbcb","bcdc","bcbd","edcb","bedc","bced","aded","adce","acdc","aedc","aced","baed","bade","bace","badc","bacd","abed","abde","abae","abad","acbe","acbd","bacb","abce","abcb","abac","adcb","abdc","abcd","decdc","dcbed","dbcbe","bdcbd","edcbe","edcbd","bedce","edbcb","becdc","bcbed","aedce","aecdc","adced","baded","badce","bacdc","baedc","baced","abded","abaed","abade","acbed","bacbe","bacbd","abcbe","abace","abacb","adcbe","adcbd","abdce","adbcb","abadc","abcdc","abcbd","abacd","edcbed","edbcbe","ebdcbd","adecdc","baedce","baecdc","abaded","bacbed","abacbe","adcbed","adbcbe","abadce","abdcbd","abacdc","abacbd"]
A5_json_core = {to_word(to_oneline(x,5)) for x in A5_json_core_list}
A5_equiv_core = set(["1", "a", "ba", "cba", "cbab", "bab", "babd", "babe", "babde", "dcbab", "dcbacb", "cbacb", "cbacbe", "dcbabd", "cbabd", "cbabed", "babed", "babede", "cbabe", "cbac", "bac", "bace", "badce", "dcbac", "cbace", "dcba", "dcbad", "dcbacd", "cbacd", "bacd", "badcd", "baced", "cbaced", "cbad", "cbaed", "cbae", "bad", "bae", "bade", "ac", "adc", "badc", "adcd", "ace", "ad", "acd", "aced", "aed", "baed", "baede", "aede", "ae", "ade", "acde", "adcde", "adce", "acede", "b", "cb", "acb", "ab", "abd", "acbd", "acbed", "abed", "abede", "abe", "abde", "acbde", "adcbde", "adcbe", "acbede", "adcb", "badcb", "bacb", "bacbe", "badcbe", "adcbd", "acbe", "cbc", "dcbc", "cbce", "dcb", "dcbd", "dcbcd", "cbd", "cbed", "cbe", "bd", "be", "bde", "c", "bc", "bdc", "bce", "bdce", "dce", "dc", "dcd", "bdcd", "ce", "d", "cd", "bcd", "cbcd", "cbced", "bced", "bdced", "dced", "adced", "badced", "adcede", "dcede", "ced", "ed", "bed", "bede", "ede", "e", "de", "cde", "bcde", "cbcde", "cbde", "dcbde", "dcbe", "dcbce", "cbede", "dcbcde", "cbcede", "bdcde", "bcede", "bdcede", "dcde", "cede"])
A5_json_reverse = {to_word(to_oneline(x[::-1],5)) for x in A5_json_core}
print(A5_json_reverse==A5_equiv_core)

# A4_list = [[0,1,2,3],[0,1,3,2],[0,2,1,3],[0,3,1,2],[1,0,2,3],[1,0,3,2],[1,2,0,3]]
# A4_word_list = [to_word(x) for x in A4_list]
# print(A4_word_list)
# print(set(A4_word_list)==set(L))