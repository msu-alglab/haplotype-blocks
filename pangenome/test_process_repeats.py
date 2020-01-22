import process_repeats as pr
import filecmp


def test_1000_19():

    # read in k and min repeat length, which will help us build the filenames
    k = "1000"             # k from de bruijn graph
    min_length = "19"    # minimum repeat length found by mummer

    # expected filenames
    filename = "yeast10.k" + k + ".fa"
    repeats_filename = "repeats.k" + k + "." + min_length + ".txt"
    pathlocs_filename = "pathlocs.k" + k + ".txt"

    # look at file and figure out how long the snp encoding and termination
    # character encodings are.
    snp_length, termination_length = pr.get_params(filename)

    f = open(filename, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    # long string of paths through cbdg
    long_string = ''.join(lines).strip()

    locs = pr.get_path_locs(termination_length, pathlocs_filename)

    # look at the repeats file and find sets of repeats
    repeats = pr.get_repeats(repeats_filename, filename)

    pr.process_repeats(
        repeats,
        long_string,
        locs,
        snp_length,
        k,
        min_length,
        test=True
    )
    file_compare = filecmp.cmp(
        "output.k1000.19_test.txt",
        "output.k1000.19.txt"
    )
    assert file_compare


def test_500_19():

    # read in k and min repeat length, which will help us build the filenames
    k = "500"             # k from de bruijn graph
    min_length = "19"    # minimum repeat length found by mummer

    # expected filenames
    filename = "yeast10.k" + k + ".fa"
    repeats_filename = "repeats.k" + k + "." + min_length + ".txt"
    pathlocs_filename = "pathlocs.k" + k + ".txt"

    # look at file and figure out how long the snp encoding and termination
    # character encodings are.
    snp_length, termination_length = pr.get_params(filename)

    f = open(filename, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    # long string of paths through cbdg
    long_string = ''.join(lines).strip()

    locs = pr.get_path_locs(termination_length, pathlocs_filename)

    # look at the repeats file and find sets of repeats
    repeats = pr.get_repeats(repeats_filename, filename)

    pr.process_repeats(
        repeats,
        long_string,
        locs,
        snp_length,
        k,
        min_length,
        test=True
    )
    file_compare = filecmp.cmp(
        "output.k500.19_test.txt",
        "output.k500.19.txt"
    )
    assert file_compare
