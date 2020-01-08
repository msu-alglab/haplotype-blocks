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
    global col
    global n
    global k
    global sets
    global nonwc
    n = len(seqs[0])
    k = len(seqs)

    global s
    s = seqs

    sets = get_sets(seqs, n, k)

    # for col in range(n):
    for col in [2]:
        nonwc = [0]*n
        """
        if col % 1000 == 1:
            print("Finding blocks at index {}".format(col + 1))
        """
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
    print("Call to DFS with i={} and rows={}".format(i, rows))
    branch_count = 0
    for b in [0, 1]:
        print("Explore b={} branch".format(b))
        rm_rows = rows.intersection(sets[i][abs(1 - b)])
        retain_rows = rows.difference(rm_rows)

        print("Rows to remove: {}".format(rm_rows))
        print("Rows to retain: {}".format(retain_rows))

        # update position i of nonwc
        for r in retain_rows:
            if s[r][i] == b:
                nonwc[i] += 1
        print("nonwc is {}".format(nonwc))
        pathOK = nonwc[i] > 0

        print("Processing nonwc for rm rows")
        # for each row that was removed, for each column other than the current
        # column, subtract 1 if the row at that column was a non-wildcard
        # character
        for r in rm_rows:
            print("r={}".format(r))
            for c in range(i + 1, col + 1):
                print("c={}".format(c))
                print("s[{}][{}]={}".format(r, c, s[r][c]))
                if s[r][c] != '*':
                    nonwc[c] -= 1
                    if nonwc[c] <= 0:
                        pathOK = False
        print("nowc is {}".format(nonwc))

        print("Is this path okay?", pathOK)
        right0 = retain_rows.intersection(sets[col + 1][0])
        right1 = retain_rows.intersection(sets[col + 1][1])
        print("Right-maximal: {}. Left-maximal: {}".format(right0, right1))
        if pathOK and bool(right0) and bool(right1):
            # BUG: we are over-counting
            branch_count += 1
            if len(retain_rows) > 1:
                if i == 0 or DFS(i - 1, retain_rows) == 1:
                    rows_to_print = [x + 1 for x in retain_rows]
                    print("Found block K={}, i={}, j={}".format(
                        rows_to_print, i + 1, col + 1))
                    """
                    g.write("K={}, i={}, j={}\n".format(
                        rows_to_print, i + 1, col + 1))
                    """
            else:
                print("Stopped pursuing because |remain_rows| == 1")

        for r in retain_rows:
            if s[r][i] == b:
                nonwc[i] -= 1

        for r in rm_rows:
            for c in range(i + 1, col + 1):
                if s[r][c] != '*':
                    nonwc[c] += 1
        print("going back up the tree. nonwc has been returned to {}".format(
            nonwc))

    print("Call to DFS with i={} and rows={} OVER".format(i, rows))
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
    blocks = [[1, 0, "*"],
              [1, 0, "*"],
              [0, 0, 0],
              [0, 0, 1]]
    find_blocks(blocks)

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
    """
