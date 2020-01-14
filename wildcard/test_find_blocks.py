from find_blocks import find_blocks
from brute_force_find_blocks import bf_find_blocks
import random


def test_small_block():
    seqs = [[1, 0, "*"],
            [1, 0, "*"],
            [0, 0, 0],
            [0, 0, 1]]
    find_blocks(seqs, "test.txt")
    f = open("test.txt")
    blocks = f.readlines()
    bf_blocks = bf_find_blocks(seqs)

    assert len(blocks) == len(bf_blocks)
    for i, b in enumerate(bf_blocks):
        assert b == blocks[i].strip()


def run_a_test(i):
    seqs = []
    for i in range(3):
        seq = []
        for j in range(4):
            val = random.random()
            if val < 0.3:
                seq.append(1)
            elif val < 0.8:
                seq.append(0)
            else:
                seq.append("*")
        if 1 not in seq and 0 not in seq:
            seq[1] = 0
        seqs.append(seq)
    f = open("testfile" + str(i) + ".txt", "w")
    for seq in seqs:
        f.write(str(seq) + "\n")

    find_blocks(seqs, "test.txt")
    f = open("test.txt")
    blocks = f.readlines()
    bf_blocks = bf_find_blocks(seqs)
    assert len(blocks) == len(bf_blocks)
    for i, b in enumerate(bf_blocks):
        assert b == blocks[i].strip()


def test_many():
    random.seed(0)
    for repetition in range(10):
        run_a_test(repetition)
