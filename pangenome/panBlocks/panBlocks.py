# TODO: make an actual python package to avoid ugly import like this
import sys
import networkx as nx
home = "/Users/lucywilliams"
sys.path.insert(1, home + "/haplotype-blocks/pangenome/panBlocks")
import block  # noqa


# functions used in computing all blocks
def check_start_and_end(repeat1):
    """We are adding this repeat to the graph, so check that it starts with an
    S and ends with an N or O."""
    assert repeat1[0] == "S"
    assert repeat1[-1] == "N" or repeat1[-1] == "O"


def get_first_line_length(filename):
    """Get the number of characters in the first line of the .fa file"""
    f = open(filename)
    first_line = f.readline()
    return len(first_line)


def check_repeats(repeat1, repeat2):
    """After cropping assert that there are no x's or y's and that the repeats
    are still the same."""
    if "X" in repeat1:
        print(repeat1, repeat2)
        raise AssertionError(
            "There is an x remaining in trimmed repeat")
    if "Y" in repeat1:
        print(repeat1, repeat2)
        raise AssertionError(
            "There is a y remaining in trimmed repeat")
    if "X" in repeat2:
        print(repeat1, repeat2)
        raise AssertionError(
            "There is an x remaining in trimmed repeat")
    if "Y" in repeat2:
        print(repeat1, repeat2)
        raise AssertionError(
            "There is a y remaining in trimmed repeat")
    if repeat1 != repeat2:
        print(repeat1)
        print(repeat2)
        raise AssertionError("repeats are not the same")


class PanBlocks:
    """A PanBlocks object builds up and holds Block objects.
    It needs a mummer repeats file, a mummer-encoded fasta file, a pathlocs
    file, and a paths file.

    Attributes:
    * mummer_filename
    * repeats_filename
    * pathlocs_filename
    * path_filename
    * r_method
    * frequencies: a dictionary with pairs of snps as keys and a list of four
    frequencies as values
    * k
    * min_length
    * path_info: dictionary mapping from name of path (name of sample) to the
    SNP nodes in the path and the nucleotide positions of those SNP nodes on
    the genome
    * SARS_COV2_path_info: dictionary mapping from name of path (name of
    sample) to the SNP nodes in the path and the nucleotide positions of those
    SNP nodes on the genome, for SARS_COV2 paths only
    * snps: a dictionary from snp nodes (both 1 and 0 versions) to selection
    coefficients
    * snp_length: length of snp encoding in mummer file
    * termination_length: length of termination character encoding in mummer
    file
    * mummer_file_string: a string containing the contents of the
    mummer-encoded fasta file of paths through the SNP graph
    * locs: a dictionary from start/end position to path ids
    * pathnames: a dictionary from start/end positions to path names
    * genome_length: length of genome for this organism
    * snp_index_offset: the number of SNPs, which is used for assigning an ID
    when creating a dotfile.
    """

    def __init__(self, mummer_filename, repeats_filename, y_0,
                 r_method="frequency"):
        # set filenames from input
        self.mummer_filename = mummer_filename
        self.repeats_filename = repeats_filename
        # infer pathlocs and path filenames from mummer filename
        self.pathlocs_filename = mummer_filename.split("mummer")[0] +\
            "pathlocs.txt"
        self.path_filename = mummer_filename.split("mummer")[0] +\
            "paths.txt"
        # infer k (de bruijn graph parameter) and repeat length from filenames
        # TODO: filenames are not consistent. this is broken.
        self.k = mummer_filename.split(".")[1][1:]
        # currently we assume that the minimum length of a repeat when running
        # mummer was always set to 19.
        self.min_length = 19
        # get path info from pathlocs and path files and  get genome length
        self.get_path_info()
        # set all_paths from path_info
        self.all_paths = [x[0] for x in self.path_info.values()]
        # get snp length and termination length in mummer encoding from mummer
        # file
        self.get_mummer_params()
        assert self.termination_length != 1
        # read in a long string from mummer-encoded fasta file of paths through
        # SNP graph
        self.get_long_string()
        # get two dictionaries: one from start/end position to path ids and one
        # from start/end position to path names
        self.get_path_locs()
        # set the method that we use to find r
        self.r_method = r_method
        if self.r_method == "frequency":
            self.create_freq_dict()
        self.y_0 = y_0

    def create_freq_dict(self):
        """Create a dictionary with pairs of SNPs as keys and a list of
        frequencies as values."""
        self.freqs = dict()
        for path_info in self.path_info.values():
            path = path_info[0]
            print("path is", path)
            for i1 in range(len(path)):
                snp1 = path[i1]
                snp1_id, snp1_value = snp1.split(":")
                for i2 in range(i1, len(path)):
                    snp2 = path[i2]
                    snp2_id, snp2_value = snp2.split(":")
                    # positions: 00 is 0 index, 01 is 1 index, 10 is 2 index,
                    # 11 is 3 index
                    pos_to_increment = 2*int(snp1_value) + int(snp2_value)
                    if (snp1_id, snp2_id) not in self.freqs:
                        self.freqs[snp1_id, snp2_id] = [0, 0, 0, 0]
                    self.freqs[snp1_id, snp2_id][pos_to_increment] += 1
        print("Frequencies dict is:")
        print(self.freqs)
        print("Path info:")
        for path_info in self.path_info.values():
            print(path_info[0])

    def get_blocks(self):
        """Return a list of all blocks in this PanBlocks object."""
        return self.blocks

    def get_num_blocks(self):
        return len(self.blocks)

    def get_mummer_params(self):
        """Look at .fa file to figure out the length of the snp encoding and the
        termination character encoding."""
        f = open(self.mummer_filename)
        f.readline()
        lines = []
        # assume that we'll find the end of the first snp and the first
        # termination character in the first 100 lines
        for i in range(100):
            lines.append(f.readline())

        ten_lines = ''.join([x.strip() for x in lines])
        first_o = ten_lines.find("O")
        first_n = ten_lines.find("N")
        # if both are not -1, take min
        if first_o != -1 and first_n != -1:
            end_of_first_snp = min(first_o, first_n)
        elif first_o != -1:
            end_of_first_snp = first_o
        else:
            end_of_first_snp = first_n
        assert end_of_first_snp != -1
        # snp goes from 1-end of first snp
        # assume that the first termination character is all X's
        first_x = ten_lines.find("X")
        start_of_term = first_x
        ten_lines = ten_lines[start_of_term:]
        next_s = ten_lines.find("S")
        # term goes from 0 to next_s - 1
        self.snp_length = end_of_first_snp - 1
        self.termination_length = next_s

    def get_long_string(self):
        """Read in mummer-encoded fasta file and get a long string of the
        contents."""
        f = open(self.mummer_filename)
        # get rid of header line
        f.readline()
        long_string = ""
        for line in f:
            long_string += line[:-1]
        self.mummer_file_string = long_string

    def get_path_locs(self):
        """Get a dictionary of start/end positions to path ids, and a separate
        dictionary of start/end positions to names."""
        f = open(self.pathlocs_filename)
        lines = f.readlines()
        locs = dict()
        names = dict()
        counter = 0
        index = 0
        for line in lines:
            if counter == 0:
                name = line.strip()
                pass
            elif counter == 1:
                start = int(line.strip()) - 1
            elif counter == 2:
                end = int(line.strip()) - 1
                # extend end to term length + 1
                end += self.termination_length
                locs[(start, end)] = index
                names[(start, end)] = name
                index += 1
            counter += 1
            counter = counter % 3
        self.locs = locs
        self.pathnames = names

    def write_blocks_to_file(self, filename):
        """Print out the blocks."""
        f = open(filename, "w")
        for b in self.blocks:
            b.write_to_file(f)
        f.close()

    def get_repeats(self):
        """Use connected components approach to find repeats."""

        get_first_line_length(self.repeats_filename)

        f = open(self.repeats_filename)
        f.readline()
        f.readline()
        lines = f.readlines()
        g = nx.Graph()
        weights = []
        print("Processing mummer repeats file...")
        counter = 0
        for line in lines:
            counter += 1
            if counter % 100000 == 0:
                print("Looking at line {} of mummer file".format(counter))
            # mummer repeats file lists the start of the first repeat, start of
            # second repeat, and length. starts are indexed beginning at 1.
            start1 = int(line.split()[0]) - 1
            start2 = int(line.split()[1]) - 1
            length = int(line.split()[2])
            self.check_this_match(start1, start2, length)
            og_length = length
            # crop this repeat to contain just a list of SNPs
            crop_start, crop_end = self.crop(start1, length,)
            # new start, end length
            start1 = start1 + crop_start
            start2 = start2 + crop_start
            length = length - crop_start - crop_end

            # if the cropped version actually contains SNPs, add to graph
            if crop_start > -1 and length > 0:
                repeat1, l, r = self.get_occ_from_long_string(start1, length)
                repeat2, l, r = self.get_occ_from_long_string(start2, length)

                g.add_edge(start1, start2, length=length)
                weights.append(length)
                assert length >= self.snp_length + 2
                check_repeats(repeat1, repeat2)
                check_start_and_end(repeat1)
            else:
                # assert that, if this repeat was long, it should have had snps
                if og_length > 50:
                    raise AssertionError(
                        "a long repeat not included")
        # for weight in weights:
        repeats = []
        num_weights = len(set(weights))
        print("There are {} distinct repeat lengths in the mummer repeats file"
              .format(num_weights))
        for weight in set(weights):
            subgraph = nx.Graph()
            print("----Processing length {} repeats".format(weight))
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
                            length = subgraph.\
                                get_edge_data(node1, node2)["length"]
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
        self.repeats = repeats

    def get_occ_from_long_string(self, start, length):
        """Look at long string and get occurrence."""
        occ = self.mummer_file_string[start:start + length]
        one_right = self.mummer_file_string[start - 1: start]
        one_left = self.mummer_file_string[start + length:start + length + 1]
        return occ, one_right, one_left

    def check_this_match(self, start1, start2, length):
        """Assert that this mummer match is equal and is right and left maximal.
        """
        # print("getting repeat1")
        repeat1, one_left1, one_right1 = self.\
            get_occ_from_long_string(start1, length)
        # print("getting repeat2")
        repeat2, one_left2, one_right2 = self.\
            get_occ_from_long_string(start2, length)
        # print("repeat1=", repeat1)
        # print("repeat2=", repeat2)
        assert repeat1 == repeat2
        assert one_left1 != one_left2
        assert one_right1 != one_right2

    def crop(self, start, length):
        """Given a start and length, crop it to include only snps.
        Will always start with `S`. """
        repeat, one_left, one_right = self.\
            get_occ_from_long_string(start, length)

        first_s = repeat.find("S")
        crop_start = first_s
        # no matter what, we should crop to the first s.
        repeat = repeat[crop_start:]
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

    def get_path_info(self):
        """Use the path file to get dictionaries from names to
        snps/positions. Also, create a list of all unique snps."""
        f = open(self.path_filename)
        counter = 0
        overall_max_position = 0
        self.path_info = dict()
        self.SARS_COV2_path_info = dict()
        self.snps = dict()
        self.snps_to_block_ids = dict()
        self.snp_locations = dict()
        for line in f:
            if counter == 0:
                # name of path
                name = line.strip()
            elif counter == 1:
                # snps in path
                snps = line.strip().split()
                for snp in snps:
                    self.snps[snp] = 0
                    self.snps_to_block_ids[snp] = 0
            elif counter == 2:
                # start positions of snps
                positions = [int(x) for x in line.strip().split()]
                max_pos = positions[-1]
                if max_pos > overall_max_position:
                    overall_max_position = max_pos
                self.path_info[name] = (snps, positions)
                if "SARS_COV2" in name:
                    self.SARS_COV2_path_info[name] = (snps, positions)
                for (snp, position) in zip(snps, positions):
                    self.snp_locations[snp] = position
            counter += 1
            counter = counter % 3
        self.genome_length = overall_max_position

    def get_distance(self, names, snps):
        """For a set of paths and set of snps, return the average genetic
        distance between start and end snps in the paths."""
        # name is the name of one of the paths with this block. snps is the
        # snps in the block.
        # need to look up the positions from path_info dict
        distances = []
        for name in names:
            info = self.path_info[name]
            path_snps = info[0]
            index = path_snps.index(snps[0])
            start_position = info[1][index]
            end_position = info[1][index + len(snps) - 1]
            distances.append(end_position - start_position)
        return sum(distances)/len(distances)

    def compute_blocks(self, test=False):
        """Put repeats into MPPHB form."""
        self.blocks = []

        print("There are {} repeats".format(len(self.repeats)))
        counter = 0
        for repeat in self.repeats:
            if counter % 1000 == 0:
                if len(self.repeats) > 1000:
                    print("Processed {} repeats".format(counter))
            counter += 1
            # get repeat info
            starts = repeat[0]
            length = repeat[1]
            start = starts[0]
            # print("Length of this repeat is {}".format(length))
            # print("-----Processing repeats at starts {}".format(starts))
            occ, one_left, one_right = self.\
                get_occ_from_long_string(start, length)
            # check that these repeats are all valid
            self.check_repeats_while_processing(starts, length, occ)
            # figure out which paths are involved in the repeat
            indices, names = self.get_path_indices_and_names(starts,
                                                             length, occ)
            # figure out which snps
            snps = self.decode_snps(occ)
            # if this is actually a block (at least 1 snp and at least two
            # paths):
            if len(snps) > 0 and len(set(indices)) > 1:

                # create Block differently based on r_method
                if self.r_method == "distance":
                    # figure out location on the genome of first and last snp
                    distance = self.get_distance(names, snps)
                    # create block with these snps, these path indices and
                    # names, and computed length and overall genome length
                    self.blocks.append(block.Block(snps,
                                                   list(set(indices)),
                                                   list(set(names)),
                                                   self.r_method,
                                                   self.y_0,
                                                   distance,
                                                   self.genome_length))
                elif self.r_method == "frequency":
                    # we have a dictoinary with keys [(startsnp, endsnp)]
                    # pointing to a list of frequencies for this pair
                    snp1 = snps[0].split(":")[0]
                    snp2 = snps[-1].split(":")[0]
                    distance = self.get_distance(names, snps)
                    self.blocks.append(block.Block(snps,
                                       list(set(indices)),
                                       list(set(names)),
                                       self.r_method,
                                       self.y_0,
                                       freqs=self.freqs[(snp1, snp2)],
                                       length=distance))

    def check_repeats_while_processing(self, starts, length, occ):
        """Given a set of repeats, check that they all match."""
        for start in starts[1:]:
            this_occ, this_one_left, this_one_right = self.\
                get_occ_from_long_string(start, length)
            if this_occ != occ:
                raise ValueError("Occurrences of a repeat are not equal")

    def get_path_indices_and_names(self, starts, length, occ):
        """Get the indices of paths from dictionary locs that are involved in this
        repeat"""
        indices = []
        names = []
        for start in starts:
            end = start + length - 1
            # print("start={}, end={}".format(start, end))
            matching_pathlocs = []
            for key in self.locs:
                if start >= key[0] and end <= key[1]:
                    index = self.locs[key]
                    name = self.pathnames[key]
                    matching_pathlocs.append(self.locs[key])
            # assert that exactly one path index matches to this repeat
            assert len(matching_pathlocs) == 1
            indices.append(str(index))
            names.append(name)
        return (indices, names)

    def decode_snps(self, occ):
        """Given a repeat, decode its A's and C's into the SNP ids."""
        # print("decoding repeat", occ)
        index = 1
        occ = occ[index:]
        snps = []
        while len(occ) >= self.snp_length + 1:
            snp = occ[:self.snp_length]
            occ = occ[self.snp_length:]
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
            # print("processed one snp")
            snps.append(snp_id + ":" + state_out)
        return snps

    def compute_k(self, snps):
        """For a list of SNPs, compute the number of paths including this
        list."""
        # TODO: make this more efficient
        count = 0
        snps_as_string = ' '.join(snps)
        for path in self.all_paths:
            # strip off 0/1 value of snp
            path = [x[:-2] for x in path]
            path_as_string = ' '.join(path)
            if snps_as_string in path_as_string:
                count += 1
        return count

    def compute_selection_coefficients(self):
        """Compute the selection coefficient for each block."""
        # print("Computing k values...")
        counter = 0
        for b in self.blocks:
            if counter % 1000 == 0:
                if len(self.blocks) > 1000:
                    print("Computed selection coefficients for {} blocks"
                          .format(counter))
            counter += 1
            snps = b.snps
            snps = [x[:-2] for x in snps]
            # print("--Processing block with snps {}".format(snps))
            k = self.compute_k(snps)
            b.compute_selection_coefficient(k)

    def set_selection_coefficients(self):
        """For each node, figure out what the max selection coefficient is."""
        # TODO: make this more efficient
        print("There are {} snps".format(len(self.snps)))
        counter = 0
        for snp in self.snps.keys():
            # get the max selection coeff for this node
            # print(f"Finding max selection coeff for snp {snp}")
            if counter % 1000 == 0:
                if len(self.snps) > 1000:
                    print("Set selection coefficients for {} snp nodes"
                          .format(counter))
            counter += 1
            max_s = 0
            for i, b in enumerate(self.blocks):
                if snp in b.snps:
                    if b.selection_coefficient > max_s:
                        max_s = b.selection_coefficient
                        block_id = i
            self.snps[snp] = max_s
            self.snps_to_block_ids[snp] = block_id

    def print_selection_coefficients(self):
        """Print out all selection coeffs."""
        print("Calling print_selection_coefficients")
        print(self.snps)

    def write_nodes(self, f):
        """Write the nodes of the SNP graph dotfile to file object f."""
        # write subgraph for all SNPs
        snps = list(set([int(x.split(":")[0]) for x in self.snps.keys()]))
        self.snp_index_offset = max(snps)
        for snp in snps:
            f.write("subgraph cluster_{}".format(snp) +
                    " { node [style=solid];\n")
            side_1_id = snp
            side_0_id = side_1_id + self.snp_index_offset
            side_1_selection_coeff = self.snps["{}:1".format(snp)]
            side_0_selection_coeff = self.snps["{}:0".format(snp)]
            side_1_red_val = f"{int(255 - side_1_selection_coeff * 255):0>2x}"
            side_0_red_val = f"{int(255 - side_0_selection_coeff * 255):0>2x}"
            side_1_color = '"#ff' + side_1_red_val + side_1_red_val + '"'
            side_0_color = '"#ff' + side_0_red_val + side_0_red_val + '"'
            f.write(
                '{} [label="1" style=filled'.format(side_1_id) +
                ' fillcolor={} color=black];\n'.format(side_1_color)
            )
            f.write(
                '{} [label="0" style=filled'.format(side_0_id) +
                ' fillcolor={} color=black];\n'.format(side_0_color)
            )
            f.write('label="SNP {}";\n}}\n\n'.format(snp))

    def write_edges(self, f):
        """Write the edges of the SNP graph dotfile to file object f."""
        f.write('subgraph base {\n')
        # a list of visually different colors. There may not be enough for
        # every path to be a distinct color.
        # generated using colorogical: http://vrl.cs.brown.edu/color
        colors = ["#72e5ef", "#6e3638", "#b3e61c", "#322cbf", "#7b9b47",
                  "#f90da0", "#47f0a3", "#bf012a", "#65f112", "#d25bfc",
                  "#096013", "#f7c5f1", "#2aa63a", "#c95e9f", "#c7dd91",
                  "#8a0458", "#fbbd13", "#657bec", "#f87945", "#2f5672",
                  "#e7ad79", "#059dc5", "#ac85a3", "#464a15"]
        index = 0
        for pathname, path_info in self.path_info.items():
            path = path_info[0]
            node_1 = path[0]
            snp1 = int(node_1.split(":")[0])
            side1 = int(node_1.split(":")[1])
            if side1 == 0:
                id1 = snp1 + self.snp_index_offset
            else:
                id1 = snp1
            first_node = True
            for node in path[1:]:
                node_2 = node
                snp2 = int(node_2.split(":")[0])
                side2 = int(node_2.split(":")[1])
                if side2 == 0:
                    id2 = snp2 + self.snp_index_offset
                else:
                    id2 = snp2
                # only label first node
                if first_node:
                    f.write(
                        '{} -> {} [label = "{}" color="{}"]\n'
                        .format(id1,
                                id2, pathname, colors[index]))
                else:
                    f.write(
                        '{} -> {} [color="{}"]\n'
                        .format(id1, id2, colors[index]))
                id1 = id2
                first_node = False
            index = (index + 1) % len(colors)
        f.write('}\n')

    def create_dot_file(self, filename):
        """Create a dotfile of the SNP graph, where nodes are colored by
        selection coefficient value."""
        f = open(filename, "w")
        # open graph
        f.write("digraph {\n\n")
        # write nodes and edges
        self.write_nodes(f)
        self.write_edges(f)
        # close graph
        f.write("}")
        f.close()

    # TODO: there are currently just two different ways for creating bed files:
    # one for yeast and one for covid19 data. This is because for yeast, we
    # just want the fewest possible different labels, so we use the first label
    # alphabetically. For covid19, we want SARS_COV2 pathnames only. This
    # should probably be done isn a nicer way.
    def create_bed_file_yeast(self, filename):
        """Create a bed file for all high-selection SNPs from a yeast data set
        (so prefer small number of sequences)"""
        f = open(filename, "w")
        for snp in self.snps:
            ordered_path_info = sorted(self.path_info.items())
            for path, path_info in ordered_path_info:
                if snp in path_info[0]:
                    index = path_info[0].index(snp)
                    position = path_info[1][index]
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(
                        path,
                        position,
                        position + 1,
                        self.snps[snp],
                        self.snps_to_block_ids[snp]
                    ))
                    break
        f.close()

    def create_bed_file(self, filename):
        """Create a bed file for snps of all paths"""
        f = open(filename, "w")
        for path, path_info in self.SARS_COV2_path_info.items():
            for i, snp in enumerate(path_info[0]):
                position = path_info[1][i]
                f.write("{}\t{}\t{}\t{}\t{}\n".format(
                    path,
                    position,
                    position + 1,
                    self.snps[snp],
                    self.snps_to_block_ids[snp]
                ))
        f.close()

    def create_block_dataset(self, filename):
        """Create a csv file containing a lot of data about each block"""
        f = open(filename, "w")
        f.write("num_snps,num_paths,bp_length,pop_size,y_t,r,s\n")
        for b in self.get_blocks():
            num_snps = len(b.snps)
            num_paths = len(b.paths)
            bp_length = b.length
            s = b.selection_coefficient
            pop_size = b.pop_size
            y_t = b.y_t
            r = b.r
            f.write("{},{},{},{},{},{},{}\n".format(
                num_snps,
                num_paths,
                bp_length,
                pop_size,
                y_t,
                r,
                s
                ))
        f.close()
