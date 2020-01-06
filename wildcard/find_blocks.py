def check_input(seqs):
    """Check that input is right format"""
    n = len(seqs[0])
    for seq in seqs:
        for i in seq:
            assert i in [0, 1, '*'], "Invalid character in input"

        assert len(seq) == n, "Not all sequences same length"


def find_blocks(seqs):
    """Find all maximal perfect wildcard haplotype blocks in a set of
    sequences"""
    check_input(seqs)
    # note: psuedocude is written with indices starting at 0. here, we will
    # start indices at 0, since that's how Python does it. But will add 1 when
    # we return the indices.
    # global variables
    global pathsets
    global col
    global n
    global k
    global sets
    pathsets = []
    print("pathsets=", pathsets)
    col = 1
    n = len(seqs[0])
    k = len(seqs)

    print("n=", n)
    print("k=", k)
    sets = get_sets(seqs, n, k)
    print(sets)

    for col in range(n):
        DFS(col, list(range(k)))


def get_sets(seqs, n, k):
    """Find the set0, set1, set* for each column."""
    sets = dict()

    for col in range(n):
        set0 = []
        set1 = []
        setw = []
        for row in range(k):
            if seqs[row][col] == 0:
                set0.append(row)
            elif seqs[row][col] == 1:
                set1.append(row)
            else:
                setw.append(row)
        this_col_sets = dict()
        this_col_sets[0] = set0
        this_col_sets[1] = set1
        this_col_sets['*'] = setw
        sets[col] = this_col_sets

    return sets


def DFS(column, rows):
    """Explore one layer deeper in the binary trie for MPWHBs at a certain
    column."""
    branch_count = 0
    for b in [0, 1]:
        pathsets.append(sets[column][b])
        print("pathsets=", pathsets)


if __name__ == "__main__":
    blocks = [[1, 0, 1, 0],
              ['*', 1, 1, 1],
              [0, '*', 1, 1]]

    find_blocks(blocks)
