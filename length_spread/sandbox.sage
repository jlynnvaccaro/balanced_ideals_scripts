# cartan = (["B",5],["A",3])
# G = WeylGroup(cartan)
# T = G.cartan_type()
# V = G.simple_reflections()
#print(V)

def switcheroo(s,start=0,end=0):
    if start == 0:
        return s[0:start] + s[end-1::-1] + s[end:]
    else:
        return s[0:start] + s[end-1:start-1:-1] + s[end:]

# def alphabet_transform(W):
#     "Returns an alphabet transformed to translate from the enumerate to the sage convention"
#     cartan = W.cartan_type().component_types()
#     alphabet = "abcdefghijklmnopqrstuvwxyz"
#     start = 0
#     end = 0
#     for t,n in cartan:
#         end += n
#         if t=="B":
#             alphabet=switcheroo(alphabet,start,end)
#         start += n
#     return alphabet

def alphabet_transform(W):
    "Returns an alphabet transformed to translate from the enumerate to the sage convention"
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    start = 0
    end = 0
    type = W.cartan_type()
    try:
        t,n = type
        end += n
        print(t,n)
        if t == "B":
            alphabet=switcheroo(alphabet,start,end)
        return alphabet
    except:
        cartan = type.component_types()
        for t,n in cartan:
            end += n
            if t=="B":
                alphabet=switcheroo(alphabet,start,end)
            start += n
        return alphabet

W = WeylGroup((["B",3]), "s")
# print(alphabet_transform(W))
# print(switcheroo("abcde",0,3))
# print(switcheroo("abcde",1,4))
# print(switcheroo("abcde",2,5))
# X = set(W.simple_reflections())
# print(next(iter(X)))
# print(cartan)

cm = CartanMatrix([[2,-1,0],[-1,2,-1],[0,-2,2]])
W = WeylGroup(cm)
for g in W.gens():
    print(g, "\n")

print("XXXXXXX")

cm2 = CartanMatrix([[2,-2,0],[-1,2,-1],[0,-1,2]])
W2 = WeylGroup(cm2)
for g in W2.gens():
    print(g, "\n")

print(W.character_table())
print(W2.character_table())
