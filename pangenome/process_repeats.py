import time
import sys
from hapBlocks.panBlocks import PanBlocks


if __name__ == "__main__":

    # start timer
    start_time = time.time()

    # get mummer filename and repeats filename from input
    mummer_filename = sys.argv[1]
    repeats_filename = sys.argv[2]
    k = mummer_filename.split(".")[1][1:]
    min_length = repeats_filename.split(".")[2]

    # create PanBlocks object and compute blocks
    b = PanBlocks(mummer_filename, repeats_filename)
    b.get_repeats()
    b.compute_blocks()
    b.write_blocks_to_file("output_file.txt")

    print("Found {} blocks".format(b.get_num_blocks()))
    # b.plot_snps_vs_paths_scatter()

    print("--- %s seconds ---" % (time.time() - start_time))

    # compute selection coeffs and make picture
    b.compute_selection_coefficients()
    b.set_selection_coefficients()
    # b.decorate_snp_graph("covid19/dotfiles/pangraph.dot")

    # generate histograms
    # b.generate_hists()

    # create dot file
    b.create_dot_file("snp_graph.dot")
