import itertools


def bf_find_blocks(seqs):
    """Find blocks with ascending j order"""

    global s
    s = seqs

    blocks = []

    for j in range(len(s[0])):
        print("Processing j={}".format(j))
        for length in range(j + 1, 0, -1):
            print("length={}".format(length))
            for prod in itertools.product("01", repeat=length):
                prod = [int(x) for x in prod]
                block = get_block(prod, length, j)
                if block:
                    row_str = [x + 1 for x in block[2]]
                    i_str = block[0] + 1
                    j_str = block[1] + 1
                    block_string = "K={}, i={}, j={}".format(
                        row_str,
                        i_str,
                        j_str)
                    blocks.append(block_string)
        print()

    print(blocks)


def get_block(prod, length, j):
    """Given a permutation, check whether it is a block with j"""
    # check which rows match
    i = j - length + 1
    matching_rows = set(range(len(s)))
    for index in range(length):
        for x, seq in enumerate(s):
            if not (seq[i + index] == prod[index] or seq[i + index] == "*"):
                if x in matching_rows:
                    matching_rows.remove(x)

    support = [False]*length
    for index in range(length):
        for x, seq in enumerate(s):
            if x in matching_rows:
                if seq[i + index] == prod[index]:
                    support[index] = True

    # check whether right-maximal
    rights = []
    rm = False
    if j + 1 == len(s[0]):
        rm = True
    else:
        for row in matching_rows:
            rights.append(s[row][j + 1])
        if 1 in rights and 0 in rights:
            rm = True
        if 1 not in rights and 0 not in rights:
            rm = True

    # check whether left-maximal
    lefts = []
    lm = False
    if i == 0:
        lm = True
    else:
        for row in matching_rows:
            lefts.append(s[row][i - 1])
        if 1 in lefts and 0 in lefts:
            lm = True
        if 1 not in lefts and 0 not in lefts:
            lm = True

    enough_rows = len(matching_rows) > 1
    all_cols_supp = False not in support

    if enough_rows and all_cols_supp and lm and rm:
        return(i, j, matching_rows)


if __name__ == "__main__":
    seqs = [[1, 0, "*"],
            [1, 0, "*"],
            [0, 0, 0],
            [0, 0, 1]]

    print("Seqs is {}".format(seqs))

    bf_find_blocks(seqs)
