import random
import sys
sys.setrecursionlimit(11000)


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
    col = 1
    n = len(seqs[0])
    k = len(seqs)

    sets = get_sets(seqs, n, k)

    for col in range(n):
        if col % 1000 == 1:
            print("Finding blocks at index {}".format(col + 1))
        DFS(col, set(range(k)))


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
        this_col_sets[0] = set(set0)
        this_col_sets[1] = set(set1)
        this_col_sets['*'] = set(setw)
        sets[col] = this_col_sets
    sets[n] = dict()
    for b in [0, 1]:
        sets[n][b] = list(range(k))

    return sets


def DFS(i, rows):
    """Explore one layer deeper in the binary trie for MPWHBs at a certain
    i."""
    branch_count = 0
    for b in [0, 1]:
        pathsets.append(sets[i][b])
        pathOK = True
        # print("Checking the {} branch".format(b))
        # print("Pathsets is currently", pathsets)
        for s in pathsets:
            if not bool(rows.intersection(s)):
                pathOK = False
                break
        rows_b = rows.intersection(sets[i][b].union(sets[i]["*"]))
        # print("Rows supporting this branch: {}".format(rows_b))
        for s in pathsets:
            if not bool(rows_b.intersection(s)):
                pathOK = False
                break
        # print("Is this path okay?", pathOK)
        right0 = rows_b.intersection(sets[col + 1][0])
        right1 = rows_b.intersection(sets[col + 1][1])
        if pathOK and bool(right0) and bool(right1) and len(rows_b) > 1:
            branch_count += 1
            if i == 0 or DFS(i - 1, rows_b) != 1:
                rows_to_print = [x + 1 for x in rows_b]
                g.write("K={}, i={}, j={}\n".format(
                    rows_to_print, i + 1, col + 1))
        pathsets.pop()
    return branch_count


def create_data(filename, k, n, wc_prop):
    """Read in binary matrix file."""
    print("Reading in file...\n")
    f = open(filename)
    seqs = []
    for line in f.readlines()[:k]:
        seq = []
        for character in line[:n]:
            if random.random() > wc_prop:
                seq.append(int(character))
            else:
                seq.append('*')
        seqs.append(seq)
    return seqs


if __name__ == "__main__":
    """
    blocks = [[1, 0, '*'],
              [1, 0, '*'],
              [0, 0, 1],
              [0, 0, '*']]
    """
    random.seed(1)
    global g

    for wc_prop in [0.0, 0.05]:
        output_name = "blocks" + str(wc_prop) + ".txt"
        g = open(output_name, "w")
        blocks = create_data(
            "/home/lucy/vcfs/output.txt",
            1000,
            10000,
            wc_prop
        )
        find_blocks(blocks)
        g.close()
