from smooth import *

class TypeA(SmoothnessChecker):
    def __init__(self,n):
        super().__init__("A{}".format(n),n+1)
        self.bad_patterns = ([3,4,1,2], [4,2,3,1])

    def letter_to_oneline(self, c, L=None):
        """Given a letter e.g. 'a', updates a list by applying the permutation."""
        i = self.alphabet.index(c)
        if L is None:
            L = [x for x in range(self.n)]
        L[i],L[i+1] = L[i+1],L[i]
        return L

# A3 principal balanced ideals
A3 = TypeA(3)
A3.schubert_smooth("caba")
A3.schubert_smooth("abcb")

# A5 principal balanced ideals
A5 = TypeA(5)
A5.schubert_smooth("abacdebcdbcb")
A5.schubert_smooth("deabcdabcaba")
A5.schubert_smooth("aedbcabaded")