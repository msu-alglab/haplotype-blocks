import time
import networkx as nx
import sys


def get_distinct_repeats(repeats_filename, repeats, long_string):
    """Find distinct repeats from repeat file."""

    # build a dict by length
    f = open(repeats_filename)
    repeat_dict = {}
    for repeat in repeats:
        start1, start2, length = repeat
        if length not in repeat_dict:
            repeat_dict[length] = []
        repeat_dict[length].append(start1)
        repeat_dict[length].append(start2)
    f.close()

    # if a length has more than one repeat, disambiguate
    clean_repeats = []
    for length in repeat_dict:
        values = repeat_dict[length]
        num_values = len(values)
        # need to disambiguate
        if num_values == 2:
            # only one repeat of this length
            values = [values[0]-1, values[1]-1]
            clean_repeats.append([values, length])
        else:
            # print("Length {} has {} elements".format(length, num_values))
            # print(values)
            occs = {}
            for s in values:
                start = s - 1
                end = s + length - 1
                occ = long_string[start:end]
                if occ not in occs:
                    occs[occ] = [start]
                else:
                    occs[occ].append(start)
            for starts in occs.values():
                clean_repeats.append([list(set(starts)), length])

    return clean_repeats


def get_occ_from_file(start, length, filename, first_line_length):
    """Look into yeast10.k.fa file for this occurrence and the character
    directly to the left and directly to the right"""
    f = open(filename)
    num_linebreaks = start // 80
    start_how_far_into_a_line = start % 80
    repeat_start = first_line_length + num_linebreaks + start
    linebreaks_in_repeat = (start_how_far_into_a_line + length) // 80
    f.seek(repeat_start)
    part_of_file = f.read(length + linebreaks_in_repeat)
    occ = "".join(part_of_file.split("\n"))
    # check whether this character before this is a linebreak
    if start % 80 == 0:
        # if it is, go back two
        f.seek(repeat_start - 2)
    else:
        # otherwise, just go back one
        f.seek(repeat_start - 1)
    one_lefts = f.read(1)
    f.seek(repeat_start + length + linebreaks_in_repeat)
    one_rights = f.read(1)
    return occ, one_lefts, one_rights


def get_first_line_length(filename):
    """Get the number of characters in the first line of the .fa file"""
    f = open(filename)
    first_line = f.readline()
    return len(first_line)


def process_repeats(repeats, locs, snp_length, k, m, filename,
                    test=False):
    if test:
        file_to_write = "output.k" + k + "." + m + "_test.txt"
    else:
        file_to_write = "output.k" + k + "." + m + ".txt"
    f = open(file_to_write, "w")
    first_line_length = get_first_line_length(filename)
    for repeat in repeats:
        starts = repeat[0]
        length = repeat[1]
        print("Length of this repeat is {}".format(length))
        start = starts[0]
        occ, one_left, one_right = get_occ_from_file(
            start,
            length,
            filename,
            first_line_length
        )
        one_lefts = set(one_left)
        one_rights = set(one_right)
        print("-----Processing repeats at starts {}".format(starts))
        for start in starts[1:]:
            # if there are three or more occurrences, we only need one mismatch
            # on right and left.
            this_occ, this_one_left, this_one_right = get_occ_from_file(
                start,
                length,
                filename,
                first_line_length
            )

            one_lefts.add(this_one_left)
            one_rights.add(this_one_right)
            if this_occ != occ:
                raise ValueError("Occurrences of a repeat are not equal")
        # since we crop the repeats to include snp characters only, these
        # repeats need not be left and right maximal.
        # if len(one_lefts) == 1:
        #    raise ValueError("Repeat is not left-maximal")
        # if len(one_rights) == 1:
        #    raise ValueError("Repeat is not right-maximal")

        # figure out which paths
        indices = []
        for start in starts:
            end = start + length - 1
            start_index = occ.find('S')
            start = start + start_index
            for key in locs:
                if start >= key[0] and end <= key[1]:
                    name = locs[key]
            indices.append(str(name))
        f.write(" ".join(indices) + "\n")
        print("Paths are:")
        print(" ".join(indices))

        # figure out which snps
        start_index = occ.find('S')
        if start_index == -1:
            f.write("No SNPS\n")
        else:
            index = start_index + 1
            occ = occ[index:]
            snps = []
            if len(occ) < snp_length + 1:
                f.write("No SNPS\n")
            while len(occ) >= snp_length + 1:
                snp = occ[:snp_length]
                occ = occ[snp_length:]
                snp = snp.replace("C", "0")
                snp = snp.replace("G", "1")
                snp_id = str(int(snp, 2))
                state = occ[0]
                occ = occ[2:]
                if state == "O":
                    state_out = "0"
                elif state == "N":
                    state_out = "1"
                else:
                    raise ValueError("State was not O or N")
                snps.append(snp_id + ":" + state_out)
            if len(set(indices)) < len(indices):
                print("Duplicate indices")
            if len(snps) > 0 and len(set(indices)) > 1:
                f.write(' '.join(snps) + "\n")
    f.close()


def get_path_locs(TERMINATION_LENGTH, pathlocs_filename):
    """Get a dictionary of start/end positions to path ids"""
    f = open(pathlocs_filename)
    lines = f.readlines()
    locs = dict()
    counter = 0
    index = 0
    for line in lines:
        if counter == 0:
            # name = line.strip()
            pass
        elif counter == 1:
            start = int(line.strip()) - 1
        elif counter == 2:
            end = int(line.strip()) - 1
            # extend end to term length + 1
            end += TERMINATION_LENGTH
            locs[(start, end)] = index
            index += 1
        counter += 1
        counter = counter % 3
    return locs


def crop(start, length, filename, first_line_length):
    """Given a start and length, crop it to include only snps."""

    repeat, one_left, one_right = get_occ_from_file(
        start,
        length,
        filename,
        first_line_length
    )

    print("Repeat looking to crop is", repeat)
    first_s = repeat.find("S")
    crop_start = first_s
    if first_s == -1:
        print("No S!")
    # no matter what, we should crop to the first s.
    repeat = repeat[crop_start:]
    print("After deleting front, repeat is", repeat)
    # find last o or n
    first_o = repeat.find("O")
    first_n = repeat.find("N")
    if first_o + first_n == -2:
        crop_end = len(repeat)
    else:
        last_o_or_n = max([i for (i, val) in enumerate(repeat) if
                           val in ("N", "O")])
        crop_end = len(repeat) - last_o_or_n - 1
    return (crop_start, crop_end)


def get_repeats(repeats_filename, filename):
    """Use connected components approach to find repeats."""

    first_line_length = get_first_line_length(filename)

    # for debugging, get long string
    f = open(filename, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    long_string = ''.join(lines).strip()

    f = open(repeats_filename)
    f.readline()
    f.readline()
    lines = f.readlines()
    g = nx.Graph()
    weights = []
    print("About to find repeats")
    counter = 0
    for line in lines:
        print("Looking at line", counter)
        counter += 1
        start1 = int(line.split()[0]) - 1
        start2 = int(line.split()[1]) - 1
        length = int(line.split()[2])
        og_length = length
        print("Looking for crops for start", start1)
        crop_start, crop_end = crop(
            start1,
            length,
            filename,
            first_line_length
        )
        start1 = start1 + crop_start
        start2 = start2 + crop_start
        length = length - crop_start - crop_end
        if length == -1:
            print("Length == -1!")
        print("After cropping")
        print("start=", start1)
        print("start=", start2)
        print("length=", length)
        print("Repeat=", long_string[start1:start1+length])
        if crop_start > -1 and length > 0:
            g.add_edge(start1, start2, length=length)
            weights.append(length)
            if "X" in long_string[start1:start1+length]:
                print(long_string[start1:start1+length])
                raise AssertionError(
                    "There is an x remaining in trimmed repeat")
            if "Y" in long_string[start1:start1+length]:
                print(long_string[start1:start1+length])
                raise AssertionError(
                    "There is a y remaining in trimmed repeat")
            if "X" in long_string[start2:start2+length]:
                print(long_string[start2:start2+length])
                raise AssertionError(
                    "There is an x remaining in trimmed repeat")
            if "Y" in long_string[start2:start2+length]:
                print(long_string[start2:start2+length])
                raise AssertionError(
                    "There is a y remaining in trimmed repeat")
        else:
            print("length {} was not a valid repeat".format(og_length))
            if og_length > 50:
                raise AssertionError(
                    "a long repeat not included")
    # for weight in weights:
    repeats = []
    for weight in set(weights):
        print("building subgraph for weight {}".format(weight))
        subgraph = nx.Graph()
        # build subgraph with ony edges this weight or greater
        for edge in g.edges(data=True):
            node1 = edge[0]
            node2 = edge[1]
            length = edge[2]["length"]
            if length >= weight:
                subgraph.add_edge(node1, node2, length=length)
        # for conn component in subgraph:
        connected_components = nx.connected_components(subgraph)
        for c in connected_components:
            # if one of the edges in the component == weight
            add_this_component = False
            for node1 in c:
                for node2 in c:
                    if subgraph.has_edge(node1, node2):
                        length = subgraph.get_edge_data(node1, node2)["length"]
                        if length == weight:
                            add_this_component = True
            if add_this_component:
                # add all edges with node1 < node2
                repeat_list = []
                for node1 in c:
                    for node2 in c:
                        if node1 > node2:
                            repeat_list.append(node1)
                            repeat_list.append(node2)
                repeats.append([list(set(repeat_list)), weight])
    return repeats


def get_params(filename):
    """Look at .fa file to figure out the length of the snp encoding and the
    termination character encoding."""
    f = open(filename)
    f.readline()
    lines = []
    # assume that we'll find the end of the first snp and the first termination
    # character in the first 10 lines (very generous assumption!)
    for i in range(10):
        lines.append(f.readline())

    ten_lines = ''.join([x.strip() for x in lines])
    first_o = ten_lines.find("O")
    first_n = ten_lines.find("N")
    end_of_first_snp = min(first_o, first_n)
    # snp goes from 1-end of first snp
    first_x = ten_lines.find("X")
    first_y = ten_lines.find("Y")
    start_of_term = min(first_x, first_y)
    ten_lines = ten_lines[start_of_term:]
    next_s = ten_lines.find("S")
    # term goes from 0 to next_s - 1
    return (end_of_first_snp - 1, next_s)


if __name__ == "__main__":
    # start timer
    start_time = time.time()

    # read in k and min repeat length, which will help us build the filenames
    k = sys.argv[1]             # k from de bruijn graph
    min_length = sys.argv[2]    # minimum repeat length found by mummer

    # expected filenames
    filename = "yeast10.k" + k + ".fa"
    repeats_filename = "repeats.k" + k + "." + min_length + ".txt"
    pathlocs_filename = "pathlocs.k" + k + ".txt"

    # look at file and figure out how long the snp encoding and termination
    # character encodings are.
    snp_length, termination_length = get_params(filename)

    # locs is dictionary with start/end indices as keys and path index as
    # values
    locs = get_path_locs(termination_length, pathlocs_filename)

    # look at the repeats file and find sets of repeats
    repeats = get_repeats(repeats_filename, filename)

    process_repeats(repeats, locs, snp_length, k, min_length,
                    filename)
    print("--- %s seconds ---" % (time.time() - start_time))
