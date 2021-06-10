class SmoothnessChecker():
    """Superclass for smoothness checkers by cartan type. Defines shared functions."""
    def __init__(self, type, n):
        """bad_patterns should be defined in the subclass."""
        self.alphabet="abcdefghijklmnopqrstuvwxyz"
        self.n = n
        self.bad_patterns = ()
        print("\n{} smoothness checker".format(type))
        print("----------------------")

    def letter_to_oneline(self, c, L=None):
        """Implement in the subclasses"""
        return NotImplemented

    def word_to_oneline(self, word):
        """Given a word e.g. 'abc', returns a permutation of length n representing the letter."""
        L = [x for x in range(self.n)]
        for c in word[::-1]:
            L = self.letter_to_oneline(c,L=L)
        return L

    def piece_to_pattern(self, L):
        """Given a 4-elt list, e.g. [1,5,2,8] converts that list to a pattern like [1,3,2,4]"""
        pattern = [0,0,0,0]
        for i in range(4):
            index = L.index(max(L))
            pattern[index] = 4-i
            L[index] = -1
        return pattern

    def schubert_smooth(self, word):
        """Checks if a schubert cell is smooth using Lakshmibai-Sandhya pattern avoidance
        Specifically, checks whether the oneline representation of the word has any bad patterns."""
        print("Checking ideal <{}>...".format(word))
        oneline = self.word_to_oneline(word)
        for a in range(self.n-3):
            for b in range(a+1, self.n-2):
                for c in range(b+1, self.n-1):
                    for d in range(c+1, self.n):
                        piece = [oneline[a], oneline[b], oneline[c], oneline[d]]
                        pattern = self.piece_to_pattern(piece)
                        if pattern in self.bad_patterns:
                            print(" * Found a bad pattern",pattern,"at indices {} {} {} {} in {}".format(a,b,c,d,oneline))
                            print(" * Balanced ideal <{}> is not smooth.".format(word))
                            return False
                        # else:
                        #     print("Good pattern:",pattern)
        print("   Balanced ideal <{}> is smooth!".format(word))
        return True