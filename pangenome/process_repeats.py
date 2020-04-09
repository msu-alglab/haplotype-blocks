import time
import sys
from hapBlocks.blocks import PanBlocks


if __name__ == "__main__":

    # start timer
    start_time = time.time()

    # get mummer filename and repeats filename from input
    mummer_filename = sys.argv[1]
    repeats_filename = sys.argv[2]
    k = mummer_filename.split(".")[1][1:]
    min_length = repeats_filename.split(".")[2]

    b = PanBlocks(mummer_filename, repeats_filename)
    b.get_repeats()
    b.compute_blocks()
    b.write_blocks_to_file("output_file.txt")
    print("--- %s seconds ---" % (time.time() - start_time))
