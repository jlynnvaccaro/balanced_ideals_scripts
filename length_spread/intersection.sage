import sys
import utils

def switcheroo(s,start=0,end=0):
    "Returns a string with the interval start:end-1 reversed"
    if start == 0:
        return s[0:start] + s[end-1::-1] + s[end:]
    else:
        return s[0:start] + s[end-1:start-1:-1] + s[end:]

def alphabet_transform(W):
    "Returns an alphabet transformed to translate from the enumerate to the sage convention"
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    start = 0
    end = 0
    type = W.cartan_type()
    try:
        t,n = type
        end += n
        if t == "B" or t == "D":
            alphabet=switcheroo(alphabet,start,end)
        return alphabet
    except:
        cartan = type.component_types()
        for t,n in cartan:
            end += n
            if t=="B" or t == "D":
                alphabet=switcheroo(alphabet,start,end)
            start += n
        return alphabet

def ab_to_s1s2(W, s):
    """Converts a Dumas-stecker string (e.g. ab) to a sage Weyl group element (e.g. s1*s2)"""
    alphabet = alphabet_transform(W)
    gens = W.simple_reflections()
    word = W.one()
    if s == "1":
        return word
    for c in s:
        word *= gens[alphabet.index(c)+1]
    return word


def s1s2_to_ab(x):
    """Converts a Weyl group element (e.g. s1*s2) to a Dumas-Stecker string ('ab')"""
    # TODO: Assumes compatible ordering
    W = x.parent()
    alphabet = alphabet_transform(W)
    elt=""
    if x==W.one():
        return "1"
    x = x.__repr__().split("*")
    for g in x:
        n = int(g[1:])
        elt += alphabet[n-1]
    return elt

def principal_ideal(x):
    """Dumas - principal ideal in W generated by x"""
    W = x.parent()
    return W.bruhat_interval(W.one(),x)

def killing_inner_product(elt1, elt2):
    """Returns the Killing inner product"""
    M1 = elt1.matrix().adjoint()
    M2 = elt2.matrix()
    # print(elt1, "\n", M1,"\n")
    # print(elt2, "\n", M2,"\n")
    print(elt1, elt2, elt1*elt2)
    print(M1*M2)
    print((M1*M2).trace())
    return (M1*M2).trace()

def cartan_matrix_entry(row_elt, column_elt):
    """Returns the row_elt,column_elt entry in the cartan matrix"""
    if row_elt == column_elt:
        return "2"
    else:
        return str(2*killing_inner_product(row_elt, column_elt) / killing_inner_product(column_elt, column_elt))

def print_abc_cartan_matrix(W):
    alphabet = alphabet_transform(W)
    num_gens = len(W.simple_reflections())
    matrix = ""
    for r in alphabet[0:num_gens]:
        gen_r = ab_to_s1s2(W,r)
        matrix += "[ "
        for c in alphabet[0:num_gens]:
            gen_c = ab_to_s1s2(W,c)
            matrix += cartan_matrix_entry(gen_r,gen_c) + " "
        matrix += "]\n"
    print(matrix)


def union_ideal(I,x):
    """Take the union of an ideal and the principal ideal generated by x."""
    L = principal_ideal(x)
    I.update(L)
    return I

def generating_set(I):
    """Given an ideal I in a Weyl group W, return its minimal generating set."""
    x = next(iter(I))
    W = x.parent()
    gens = W.simple_reflections()
    I_gens = set()
    for x in I:
        is_gen = True
        for g in gens:
            # If g*x is an upper cover of x, then it's not a generator.
            if g*x in I and (g*x).length() > x.length():
                is_gen = False
                break
            if x*g in I and (x*g).length() > x.length():
                is_gen = False
                break
        if is_gen:
            I_gens.add(x)
    return I_gens
        

def perp(I):
    """Given a Weyl group W and an ideal I, calculates w0*I_C (longest elt times the complement of I)"""
    x = next(iter(I))
    W = x.parent()
    w0 = W.long_element()
    I_C = set(principal_ideal(w0)) - I
    return {w0*x for x in I_C}


def print_core_hull(infile):
    d = utils.load_json_data(infile)
    summands = d["summands"]
    cartan = []
    for s in summands:
        print(s)
        L = [s[0],int(s[1:])]
        cartan.append(L)

    G = WeylGroup(cartan, "s")# * WeylGroup(cartan, "s")


    core = set(principal_ideal(G.long_element()))
    hull = set()
    Ideals = d["ideals"]

    for I in Ideals:
        b_set = set()
        for g in I:
            word = ab_to_s1s2(G,g)
            b_set = union_ideal(b_set, word)
        core = core.intersection(b_set)
        hull = hull.union(b_set)


    
    print("Size of the core:",len(core))
    print("Size of the hull:",len(hull))
    abc_core = {s1s2_to_ab(sword) for sword in core}
    abc_hull = {s1s2_to_ab(sword) for sword in hull}

    print("Core elements:",abc_core)
    print("Hull elements:",abc_hull)
    abc_coregen = {s1s2_to_ab(sword) for sword in generating_set(core)}
    abc_hullgen = {s1s2_to_ab(sword) for sword in generating_set(hull)}

    print("Core generators:",abc_coregen)
    print("Hull generators:",abc_hullgen)
    core_perp = perp(core)
    hull_perp = perp(hull)
    
    print("Does core-perp=hull?",core_perp == hull)
    print("Does hull-perp=core?",hull_perp == core)

# Actual script starts here...
# TODO: maybe move the functions to utils?

if __name__=="__main__":

    # if len(sys.argv) < 2:
    #     raise TypeError("Requires an input arg, e.g. 'sage intersection.sage A3.json'")

    # infile = sys.argv[1]
    # print_core_hull(infile)
    G = WeylGroup("A3","s")
    print_abc_cartan_matrix(G)
