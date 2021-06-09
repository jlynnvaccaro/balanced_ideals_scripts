global alphabet
global bad_patterns
alphabet="abcdefghijklmnopqrstuvwxyz"
bad_patterns = [[3,4,1,2], [4,2,3,1]]

def letter_to_oneline(c, L=None, n=0):
    """Given a letter e.g. 'a', updates a list by applying the permutation."""
    i = alphabet.index(c)
    if L is None:
        L = [x for x in range(n)]
    L[i],L[i+1] = L[i+1],L[i]
    return L

def word_to_oneline(word, n):
    """Given a word e.g. 'abc', returns a permutation of length n representing the letter."""
    L = [x for x in range(n)]
    for c in word[::-1]:
        L = letter_to_oneline(c,L=L)
    return L

def piece_to_pattern(L):
    """Given a 4-elt list, e.g. [1,5,2,8] converts that list to a pattern like [1,3,2,4]"""
    pattern = [0,0,0,0]
    for i in range(4):
        index = L.index(max(L))
        pattern[index] = 4-i
        L[index] = -1
    return pattern

def schubert_smooth(word,n=10):
    """Checks if a schubert cell is smooth using Lakshmibai-Sandhya pattern avoidance
    Specifically, checks whether the oneline representation of the word has any bad patterns."""
    print("Checking ideal <{}>...".format(word))
    oneline = word_to_oneline(word,n)
    for a in range(n-3):
        for b in range(a+1, n-2):
            for c in range(b+1, n-1):
                for d in range(c+1, n):
                    piece = [oneline[a], oneline[b], oneline[c], oneline[d]]
                    pattern = piece_to_pattern(piece)
                    if pattern in bad_patterns:
                        print(" * Found a bad pattern",pattern,"at indices {} {} {} {} in {}".format(a,b,c,d,oneline))
                        print(" * Balanced ideal <{}> is not smooth.\n".format(word))
                        return False
                    # else:
                    #     print("Good pattern:",pattern)
    print("   Balanced ideal <{}> is smooth!\n".format(word))
    return True

if __name__=="__main__":

    # A3 principal balanced ideals
    schubert_smooth("caba",4)
    schubert_smooth("abcb",4)

    # A5 principal balanced ideals
    schubert_smooth("abacdebcdbcb",6)
    schubert_smooth("deabcdabcaba",6)
    schubert_smooth("aedbcabaded",6)