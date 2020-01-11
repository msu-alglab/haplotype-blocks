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


def find_blocks(seqs, outputfile):
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
    global cons_cnt

    n = len(seqs[0])
    k = len(seqs)

    global s
    s = seqs

    sets = get_sets(seqs, n, k)

    global g
    g = open(outputfile, "w")

    # for col in range(n):
    for col in range(n):
        cons_cnt = [0]*n
        if col % 1000 == 1:
            print("Finding blocks at index {}".format(col + 1))
        DFS(col, set(range(k)))

    g.close()


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
    sets[n]['*'] = list(range(k))

    return sets


def DFS(i, rows):
    """Explore one layer deeper in the binary trie for MPWHBs at a certain
    i."""
    # print("Call to DFS with i={} and rows={}".format(i, rows))
    branch_count = 0
    for b in [0, 1]:
        # print("Explore b={} branch".format(b))
        kp = rows.intersection(sets[i][b].union(sets[i]['*']))
        rm = rows.intersection(sets[i][abs(1-b)])

        # print("Rows to remove: {}".format(rm))
        # print("Rows to retain: {}".format(kp))

        # update position i of nonwc
        for r in kp:
            if s[r][i] == b:
                cons_cnt[i] += 1
        # print("cons_cnt is {}".format(cons_cnt))
        consOK = cons_cnt[i] > 0

        # print("Processing nonwc for rm rows")
        # for each row that was removed, for each column other than the current
        # column, subtract 1 if the row at that column was a non-wildcard
        # character
        for r in rm:
            # print("r={}".format(r))
            for c in range(i + 1, col + 1):
                # print("c={}".format(c))
                # print("s[{}][{}]={}".format(r, c, s[r][c]))
                if s[r][c] != '*':
                    cons_cnt[c] -= 1
                    if cons_cnt[c] <= 0:
                        consOK = False
        # print("cons_cnt is {}".format(cons_cnt))

        # print("Is this path okay?", consOK)
        right0 = bool(kp.intersection(sets[col + 1][0]))
        right1 = bool(kp.intersection(sets[col + 1][1]))
        rightwc = kp.issubset(sets[col + 1]['*'])
        # print("Right0={}, right1={}, rightallwc={}".format(
        #    right0,
        #    right1,
        #    rightwc))
        r_max = col == n or right0 and right1 or rightwc
        if consOK and r_max:
            branch_count += 1
            # print("Branch count increment, so bc={}".format(branch_count))
            if len(kp) > 1:
                if i == 0 or DFS(i - 1, kp) != 1:
                    rows_to_print = [x + 1 for x in kp]
                    # print("Found block K={}, i={}, j={}".format(
                    #     rows_to_print, i + 1, col + 1))
                    g.write("K={}, i={}, j={}\n".format(
                        rows_to_print, i + 1, col + 1))
            else:
                # print("Stopped pursuing because |remain_rows| == 1")
                pass

        for r in kp:
            if s[r][i] == b:
                cons_cnt[i] -= 1

        for r in rm:
            for c in range(i + 1, col + 1):
                if s[r][c] != '*':
                    cons_cnt[c] += 1
        # print("going back up the tree. nonwc has been returned to {}".format(
        #     cons_cnt))

    # print("Call to DFS with i={} and rows={} OVER".format(i, rows))
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

    random.seed(1)

    for wc_prop in [0.05]:
        output_name = "blocks" + str(wc_prop) + ".txt"
        blocks = create_data(
            "/home/lucy/vcfs/output.txt",
            100,
            10000,
            wc_prop
        )
        find_blocks(blocks, output_name)
