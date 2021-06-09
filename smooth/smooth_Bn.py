from smooth import *

class TypeB(SmoothnessChecker):
    def __init__(self, n):
        super().__init__(n)
        self.bad_patterns = ([-2,-1],[1,2,-3],[1,-2,-3],[-1,2,-3],[2,-1,-3],[-2,1,-3],[3,-2,1],[2,-4,3,1],[-2,-4,3,1],[3,4,1,2],[3,4,-1,2],[-3,4,1,2],[4,1,3,-2],[4,-1,3,-2],[4,2,3,1],[4,2,3,-1],[-4,2,3,1])

    def letter_to_oneline(self, c, L=None):
        """Given a letter e.g. 'a', updates a list by applying the SIGNED permutation."""
        if L is None:
            L = [x for x in range(self.n)]
        if c == 'a':
            L[0] = -L[0]
        else:
            i = self.alphabet.index(c)
            L[i-1],L[i] = L[i-1],L[i]
        return L

B2 = TypeB(2)
B2.schubert_smooth("ba")
B2.schubert_smooth("ab")

B3 = TypeB(3)
B3.schubert_smooth("bcabab")
B3.schubert_smooth("cbabcb")
B3.schubert_smooth("ababcb")

B4 = TypeB(4)
B4.schubert_smooth("dcbadcbacdc")
B4.schubert_smooth("dcababcdbcb")
B4.schubert_smooth("bcdabacbabcb")
B4.schubert_smooth("abcababcdbcb")

B5 = TypeB(5)
B5.schubert_smooth("edabcababcdebcdbcb")
B5.schubert_smooth("cedcebcdabcababded")
B5.schubert_smooth("dcedbcababcdebcdbcb")
B5.schubert_smooth("bcdeabacbdcababcdbcb")
B5.schubert_smooth("abcdabcababcdebcdbcb")