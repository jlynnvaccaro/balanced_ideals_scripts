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

def to_str_list(L):
    """Takes a list of nodes and returns a list of strings without modifying the first list"""
    new_L = []
    for node in L:
        new_L.append(node.word)
    return new_L

class matrix:
    def __init__(self, data=None):
        if data is None:
            self.data = None
            return
        if isinstance(data[0],list):
            self.data = data
        elif isinstance(data[0],int):
            self.data = [[int(i==j) for i in range(len(data))] for j in data]
        self.empty = False
    def data_str(self):
        return [[str(i) for i in row] for row in self.data]

    def __repr__(self):
        str_list = [" ".join(row) for row in self.data_str()]
        return "[[" +"]\n [".join(str_list) + "]]"
    def sum(self):
        """Returns the sigma permutation matrix"""
        num_rows = len(self.data)
        copy_data = [[i for i in row] for row in self.data]
        for i in range(num_rows):
            for j in range(num_rows):
                copy_data[i][j] = sum([sum([elt for elt in row[:j+1]]) for row in self.data[:i+1]])
        return matrix(copy_data)

    def __lt__(self, other):
        """This is really <= but w.e"""
        num_rows = len(self.data)            
        truth = [all([self.data[i][j]<=other.data[i][j] for i in range(num_rows)]) for j in range(num_rows)]
        return all(truth)
    def __gt__(self,other):
        return other<self
    def total_sum(self):
        return sum([sum(row) for row in self.data])
    def __sub__(self, other):
        new_data = [[i-j for i,j in zip(row_i,row_j)] for row_i,row_j in zip(self.data,other.data)]
        return matrix(new_data)
    def __add__(self, other):
        """Returns the elt-wise sum of two matrices"""
        new_data = [[i+j for i,j in zip(row_i,row_j)] for row_i,row_j in zip(self.data,other.data)]
        return matrix(new_data)
    def dot(self):
        """Creates a diff matrix from 1 to w0"""
        num_rows = len(self.data)
        # print("w0 matrix!!! shift:",shift)
        # print(matrix([[max(i+j-shift,0) for i in range(num_rows-1)] for j in range(num_rows-1)]))
        # shift = num_rows-2
        # return matrix([[self.data[j][i] - max(i+j-shift,0) for i in range(num_rows-1)] for j in range(num_rows-1)])
        # print(matrix([[min(i+1,j+1) for i in range(num_rows-1)] for j in range(num_rows-1)]))
        return matrix([[min(i+1,j+1)-self.data[j][i] for i in range(num_rows-1)] for j in range(num_rows-1)])
    def flip(self):
        return matrix([[x for x in row[::-1]] for row in self.data])
    def dotsum(self, other):
        """Returns the elt-wise sum of two matrices"""
        return self.dot().flip() + other.dot()
    def concat_right(self, other):
        if self.data is None:
            return matrix(other.data_copy())
        new_data = [self_row + [-1] + other_row for self_row,other_row in zip(self.data,other.data)]
        return matrix(new_data)
    def data_copy(self):
        if self.data is None:
            return []
        return [[x for x in row] for row in self.data]
    def concat_down(self,other):
        if self.data is None:
            return matrix(other.data_copy())
        self_data = self.data_copy()
        other_data = other.data_copy()

        while(len(self_data[0])<len(other_data[0])):
            for i in range(len(self_data)):
                self_data[i] = [-1]+self_data[i]+[-1]
        while len(self_data[0]) > len(other_data[0])+1:
            for i in range(len(other_data)):
                other_data[i] = [-1]+other_data[i]+[-1]  
        if len(self_data[0]) > len(other_data[0]):
            for i in range(len(other_data)):
                other_data[i] = other_data[i]+[-1] 
        neg1 = [-1 for _ in self_data[0]]
        self_data += [neg1]
        for row in other_data:
            self_data += [row]
        return matrix(self_data)



class equiv_node:
    def __init__(self,oneline):
        self.oneline=oneline
        self.oneline_len = len(oneline)
        self.word = to_word(oneline)
        if self.word =="1":
            self.word_len = 0
        else:
            self.word_len = len(self.word)
        self.equiv_children = []
        self.ascend_children = []
        self.sym_children = []
        self.equiv_checked = False
        self.ascend_checked = False
        self.sym_checked = False
        self.find_triples()
    def up_set(self,n):
        self.oneline += [i for i in range(len(self.oneline),n+1)]
        self.equiv_checked = False
        self.ascend_checked = False
        self.sym_checked = False
    def set_parent(self,parent, type=""):
        parent.add_child(self,type)
        if type=="equiv":
            self.equiv_parent=parent
        if type=="ascend":
            self.ascend_parent=parent
        if type=="sym":
            self.sym_parent=parent
    def add_child(self,child, type=""):
        if type == "equiv":
            self.equiv_children.append(child)
        if type == "ascend":
            self.ascend_children.append(child)
        if type == "sym":
            self.sym_children.append(child)
    def needs_checking(self, type=""):
        if type == "equiv":
            return not self.equiv_checked
        if type == "ascend":
            return not self.ascend_checked
        if type == "sym":
            return not self.sym_checked
    def get_children(self, type=""):
        if type == "equiv":
            return self.equiv_children
        if type == "ascend":
            return self.ascend_children
        if type == "sym":
            return self.sym_children
    def __eq__(self, other):
        return self.oneline==other.oneline
    def __repr__(self):
        return self.word
    def __str__(self):
        return "{}:{}".format(self.word,self.num_triples)
    def find_triples(self):
        num_triples = 0
        for i in range(self.oneline_len-2):
            if self.oneline[i]<self.oneline[i+1]<self.oneline[i+2]:
                num_triples+=1
        self.num_triples = num_triples
    def w0(self):
        w0 = [self.oneline_len - 1 - i for i in range(self.oneline_len)]
        new_oneline = [w0[i] for i in self.oneline]
        return equiv_node(new_oneline)
    def matrix_permutation(self):
        return matrix(self.oneline)
    def matrix_sigma(self):
        I = self.matrix_permutation()
        return I.sum()


class equiv_tree:
    def __init__(self,n):
        self.root = equiv_node(to_oneline("1",n))
        self.n = n
        self.range = n-1
    def up_cartan(self,n):
        """Changes the size upwards to cartan An"""
        if n==self.n:
            return
        if n<self.n:
            return
        self.n = n
        self.range = n-1
        nodes = self.nodes()
        for node in nodes:
            node.up_set(n)

    def nodes(self,type="",L=None,elt=None):
        """Returns a list of equivalence class nodes"""
        if L is None:
            L = []
            return_L = True
        else: 
            return_L = False
        if elt is None:
            elt = self.root
        L.append(elt)
        for child in elt.get_children(type):
            self.nodes(type,L,child)
        if return_L:
            return L
    def nodes_to_check(self,type="",L=None,elt=None):
        """Returns a list of equivalence class leafs"""
        if L is None:
            L=[]
            if elt is None:
                elt=self.root
            self.nodes(type,L,elt)
        return [x for x in L if x.needs_checking(type)]
    def make_equiv_class(self):
        """Generates the nodes related to the equivalence class"""
        num_loops = 0
        nodes = self.nodes(type="equiv")
        check_nodes = self.nodes_to_check(type="equiv",L=nodes)
        while True:
            new_nodes = []
            num_loops += 1            
            if len(check_nodes) == 0:
                return num_loops
            for node in check_nodes:
                for i in range(self.range):
                    if is_increasing(node.oneline[i:i+3]) or is_increasing(new_pattern_1(node.oneline[i:i+3])):
                        new_node = equiv_node(new_pattern_1(node.oneline,i))
                        if new_node not in nodes and new_node not in new_nodes:
                            new_node.set_parent(node,"equiv")
                            nodes.append(new_node)
                            new_nodes.append(new_node)
                    if is_increasing(node.oneline[i:i+3]) or is_increasing(new_pattern_2(node.oneline[i:i+3])):
                        new_node = equiv_node(new_pattern_2(node.oneline,i))
                        if new_node not in nodes and new_node not in new_nodes:
                            new_node.set_parent(node,"equiv")
                            nodes.append(new_node)
                            new_nodes.append(new_node)
                node.checked = True
            check_nodes = new_nodes
    
    def make_ascending_equiv_class(self):
        """Generates the nodes related to the equivalence class"""
        num_loops = 0
        nodes = self.nodes(type="ascend")
        check_nodes = self.nodes_to_check(type="ascend",L=nodes)
        while True:
            new_nodes = []
            num_loops += 1            
            if len(check_nodes) == 0:
                return num_loops
            for node in check_nodes:
                for i in range(self.range):
                    if is_increasing(node.oneline[i:i+3]):
                        new_node = equiv_node(new_pattern_1(node.oneline,i))
                        if new_node not in nodes and new_node not in new_nodes:
                            new_node.set_parent(node,"ascend")
                            nodes.append(new_node)
                            new_nodes.append(new_node)
                    if is_increasing(node.oneline[i:i+3]):
                        new_node = equiv_node(new_pattern_2(node.oneline,i))
                        if new_node not in nodes and new_node not in new_nodes:
                            new_node.set_parent(node,"ascend")
                            nodes.append(new_node)
                            new_nodes.append(new_node)
                node.checked = True
            check_nodes = new_nodes

    def make_symmetric_group(self):
        """Generates the nodes related to the equivalence class"""
        num_loops = 0
        nodes = self.nodes(type="sym")
        start_num_nodes = len(nodes)
        check_nodes = self.nodes_to_check(type="sym",L=nodes)
        while True:
            new_nodes = []
            num_loops += 1            
            if len(check_nodes) == 0:
                end_len_nodes = len(nodes)
                return num_loops
            for node in check_nodes:
                for i in range(self.range):
                    new_node = equiv_node(new_pattern_1(node.oneline,i))
                    if new_node not in nodes and new_node not in new_nodes:
                        new_node.set_parent(node,"sym")
                        nodes.append(new_node)
                        new_nodes.append(new_node)    
                    new_node = equiv_node(new_pattern_2(node.oneline,i))
                    if new_node not in nodes and new_node not in new_nodes:
                        new_node.set_parent(node,"sym")
                        nodes.append(new_node)
                        new_nodes.append(new_node)
                node.checked = True
            check_nodes = new_nodes

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


if __name__=="__main__":
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