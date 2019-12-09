import scipy.special
import networkx as nx
import sys

def get_distinct_repeats(repeats_filename, repeats, long_string):
    """Find distinct repeats from repeat file."""

    # build a dict by length
    f = open(repeats_filename)
    lines = f.readlines()
    repeat_dict = {}
    for repeat in repeats:
        start1, start2, length = repeat
        if not length in repeat_dict:
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
            #print("Length {} has {} elements".format(length, num_values))
            #print(values)
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

def process_repeats(repeats, long_string, locs, snp_length, k, m):
    f = open("output.k"+k+"."+m+".txt", "w")
    for repeat in repeats:
        starts = repeat[0]
        length = repeat[1]
        #print("Length of this repeat is {}".format(length))
        start = starts[0]
        end = start + length
        occ = long_string[start:end]
        one_left = long_string[start - 1]
        one_right = long_string[end]
        one_lefts = set(one_left)
        one_rights = set(one_right)
        print("-----Processing repeats at starts {}".format(starts))
        for start in starts[1:]:
            # if there are three or more occurrences, we only need one mismatch on
            # right and left.
            end = start + length
            this_occ = long_string[start:end]
            this_one_left = long_string[start - 1]
            this_one_right = long_string[end]
            one_lefts.add(this_one_left)
            one_rights.add(this_one_right)
            if this_occ != occ:
                raise ValueError("Occurrences of a repeat are not equal")
        # since we crop the repeats to include snp characters only, these
        # repeats need not be left and right maximal.
        #if len(one_lefts) == 1:
        #    raise ValueError("Repeat is not left-maximal")
        #if len(one_rights) == 1:
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
            if len(snps) > 0:
                f.write(' '.join(snps) + "\n")
    f.close()


def get_path_locs(TERMINATION_LENGTH, pathlocs_filename):
    f = open(pathlocs_filename)
    lines = f.readlines()
    locs = dict()
    counter = 0
    index = 0
    for line in lines:
        if counter == 0:
            name = line.strip()
        elif counter == 1:
            start = int(line.strip()) - 1
        elif counter == 2:
            end = int(line.strip()) - 1
            # extend end to term length + 1
            end += TERMINATION_LENGTH
            locs[(start,end)] = index
            index += 1
        counter += 1
        counter = counter % 3
    return locs


def crop(start, length):
    """Given a start and length, crop it to include only snps."""
    f = open(filename, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    long_string = ''.join(lines).strip()

    repeat = long_string[start:start + length]
    first_s = repeat.find("S")
    first_x = repeat.find("X")
    first_y = repeat.find("Y")
    first_term_char  = min(first_x, first_y)
    contains_term_chars = first_x + first_y != -2
    if contains_term_chars and first_s >= first_term_char:
        # starts with term chars
        crop_start = first_s
    else:
        crop_start = 0
    repeat = repeat[crop_start:]
    first_x = repeat.find("X")
    first_y = repeat.find("Y")
    contains_term_chars = first_x + first_y != -2
    first_term_char = min(first_x, first_y)
    if contains_term_chars:
        crop_end = len(repeat) - first_term_char
    else:
        crop_end = 0
    return (crop_start, crop_end)


def get_repeats(repeats_filename):
    """Use connected components approach to find repeats."""
    # for debugging, get long string
    f = open("yeast10.k1000.fa", "r")
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
        crop_start, crop_end = crop(start1, length)
        start1 = start1 + crop_start
        start2 = start2 + crop_start
        length = length - crop_start - crop_end
        g.add_edge(start1, start2, length=length)
        weights.append(length)
        """
        if "X" in long_string[start1:start1+length]:
            print(long_string[start1:start1+length])
            raise AssertionError("There is an x remaining in trimmed repeat")
        if "Y" in long_string[start1:start1+length]:
            print(long_string[start1:start1+length])
            raise AssertionError("There is a y remaining in trimmed repeat")
        if "X" in long_string[start2:start2+length]:
            print(long_string[start2:start2+length])
            raise AssertionError("There is an x remaining in trimmed repeat")
        if "Y" in long_string[start2:start2+length]:
            print(long_string[start2:start2+length])
            raise AssertionError("There is a y remaining in trimmed repeat")
        """
    # for weight in weights:
    repeats = []
    for weight in set(weights):
        print("building subgraph for weight {}".format(weight))
        subgraph = nx.Graph()
        # build subgraph with ony edges this weight or greater
        for edge in g.edges(data=True):
            node1 = edge[0]
            node2 = edge[1]
            l = edge[2]["length"]
            if l >= weight:
                subgraph.add_edge(node1, node2, length=l)
        # for conn component in subgraph:
        connected_components = nx.connected_components(subgraph)
        for c in connected_components:
            # if one of the edges in the component == weight
            add_this_component = False
            for node1 in c:
                for node2 in c:
                    if subgraph.has_edge(node1, node2):
                        l = subgraph.get_edge_data(node1, node2)["length"]
                        if l == weight:
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
    lines = f.readlines()
    long_string = ''.join([x.strip() for x in lines])
    first_o = long_string.find("O")
    first_n = long_string.find("N")
    end_of_first_snp = min(first_o, first_n)
    # snp goes from 1-end of first snp
    first_x = long_string.find("X")
    first_y = long_string.find("Y")
    start_of_term = min(first_x, first_y)
    long_string = long_string[start_of_term:]
    next_s = long_string.find("S")
    # term goes from 0 to next_s - 1
    return (end_of_first_snp - 1, next_s)


if __name__ == "__main__":
    k = sys.argv[1]
    min_length = sys.argv[2]
    assert k in  ["100", "1000"]
    assert min_length == "30"

    filename = "yeast10.k" + k + ".fa"
    repeats_filename = "repeats.k" + k + "." + min_length + ".txt"
    pathlocs_filename = "pathlocs.k" + k + ".txt"

    snp_length, termination_length = get_params(filename)

    f = open(filename, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    long_string = ''.join(lines).strip()

    locs = get_path_locs(termination_length, pathlocs_filename)
    repeats = get_repeats(repeats_filename)
    #clean_repeats = get_distinct_repeats(repeats_filename, repeats, long_string)

    process_repeats(repeats, long_string, locs, snp_length, k, min_length)
