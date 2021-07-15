
def tuple_equal(T1,T2):
    return set(T1) == set(T2)

def set_tuple_equal(S,T):
    return S == set(T)

def is_slim(S):
    """Returns true if a subset of W is slim. Otherwise, returns false."""
    x = next(iter(S))
    W = x.parent()
    w0 = W.long_element()
    for x in S:
        if w0*x in S:
            return False
    return True

def principal_ideal(x):
    """Dumas - principal ideal in W generated by x"""
    W = x.parent()
    return W.bruhat_interval(W.one(),x)

def set_in_list(L,S):
    for t in L:
        if set_tuple_equal(S,t):
            return True
    return False

def slims_from_set_recursive(I, d=None,max_depth=3,depth=1):
    """Blah blah blah"""
    d[int(depth)].append(tuple(I))
    if depth>=max_depth:
        return
    upper_covers = set()
    for x in I:
        upper_covers.update(set(x.upper_covers()).difference(I))
    for u in upper_covers:
        new_I = set()
        new_I.update(I)
        new_I.update(principal_ideal(u))
        new_depth = len(new_I)
        if is_slim(new_I) and new_depth<=max_depth:
            if not set_in_list(d[new_depth],new_I):
                # print(new_I)
                slims_from_set_recursive(new_I,d,max_depth, new_depth)
    
import sys

W = WeylGroup(sys.argv[1],"s")
I_0 = {W.first()}
max_depth = len(W)/2
d = {}
for i in range(1,max_depth+1):
    d[i] = []
slims_from_set_recursive(I_0,d,max_depth=max_depth)
print(d)
print("Len: # slim ideals")
for i in range(1,max_depth+1):
    print("{}: {}".format(i,len(d[i])))
