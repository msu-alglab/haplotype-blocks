import scipy.special
import networkx as nx

def get_distinct_repeats(long_string):
    """Find distinct repeats from repeat file."""

    # build a dict by length
    f = open("repeats.txt", "r")
    lines = f.readlines()
    repeats = []
    repeat_dict = {}
    for line in lines[2:]:
        start1 = int(line.split()[0])
        start2 = int(line.split()[1])
        length = int(line.split()[2])
        repeats.append((start1,start2, length))
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

def process_repeats(repeats, long_string, locs, snp_length):
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
        if len(one_lefts) == 1:
            raise ValueError("Repeat is not left-maximal")
        if len(one_rights) == 1:
            raise ValueError("Repeat is not right-maximal")

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
        print(" ".join(indices))

        # figure out which snps
        start_index = occ.find('S')
        if start_index == -1:
            print("No SNPS")
        else:
            index = start_index + 1
            occ = occ[index:]
            snps = []
            if len(occ) < snp_length + 1:
                print("No SPS")
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
                print(' '.join(snps))


def get_path_locs(TERMINATION_LENGTH):
    f = open("pathlocs.txt")
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
            end += TERMINATION_LENGTH
            locs[(start,end)] = index
            index += 1
        counter += 1
        counter = counter % 3
    return locs


def get_repeats():
    """Use connected components approach to find repeats."""
    f = open("repeats.txt", "r")
    f.readline()
    f.readline()
    lines = f.readlines()
    g = nx.Graph()
    for line in lines:
        start1 = int(line.split()[0])
        start2 = int(line.split()[1])
        length = int(line.split()[2])
        g.add_edge(start1, start2, length=length)
    print("# connected components:",
        nx.number_connected_components(g))
    connected_components = nx.connected_components(g)
    repeats = []
    for c in connected_components:
        print("#### Processing connected component", c)
        weights = []
        for node1 in c:
            for node2 in c:
                if node1 < node2:
                    if g.has_edge(node1, node2):
                        attributes = g.get_edge_data(node1, node2)
                        weights.append(attributes["length"])
        weights = list(set(weights))
        for weight in weights:
            hap_block = []
            for node1 in c:
                for node2 in c:
                    if node1 < node2:
                        if g.has_edge(node1, node2):
                            if g.get_edge_data(node1, node2)["length"]\
                                             >= weight:
                                hap_block.append(node1)
                                hap_block.append(node2)
            hap_block = list(set(hap_block))
            print("hap block for weight {}: {}".format(
                        weight,
                        hap_block))


if __name__ == "__main__":
    SNP_LENGTH = 18
    FILENAME =  "yeast10_k1000.fa"
    TERMINATION_LENGTH = 13
    SNP_ENCODING_LENGTH = 11

    f = open(FILENAME, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    long_string = ''.join(lines).strip()

    locs = get_path_locs(TERMINATION_LENGTH)
    get_repeats()
    #clean_repeats = get_distinct_repeats(long_string)

    #process_repeats(clean_repeats, long_string, locs, SNP_ENCODING_LENGTH)
