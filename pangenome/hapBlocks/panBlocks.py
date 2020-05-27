import networkx as nx


class Block:
    """A Block object contains the SNP nodes (ID and 0/1), paths (by both ID
    and name) for a maximal pangenome haplotype block. Based on this
    information, it can compute the selection coefficient."""

    def __init__(self, snps, paths, pathnames):
        self.snps = snps
        self.paths = paths
        self.pathnames = pathnames

    def write_to_file(self, f):
        """Write this block's info to passed file."""
        f.write(" ".join(self.paths) + "\n")
        f.write(" ".join(self.snps) + "\n")

    def compute_selection_coefficient(self, k):
        """Compute the selection coefficient for this block."""
        self.selection_coefficient = len(self.snps) / k

    def compute_recombination_frequency(self, method=None):
        """Compute the recombination frequency for this block."""
        # viruses do not reproduce sexually, so I'm setting this to 0 for now.
        if method is None:
            # no recombination
            self.r = 0
        elif method == "distance":
            # compute recombination factor based on number nucleotides
            pass


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
    * k
    * min_length
    * path_info: dictionary mapping from name of path (name of sample) to the
    SNP nodes in the path and the nucleotide positions of those SNP nodes on
    the genome
    * snp_length: length of snp encoding in mummer file
    * termination_length: length of termination character encoding in mummer
    file
    * mummer_file_string: a string containing the contents of the
    mummer-encoded fasta file of paths through the SNP graph
    * locs: a dictionary from start/end position to path ids
    * pathnames: a dictionary from start/end positions to path names
    """

    def __init__(self, mummer_filename, repeats_filename):
        # set filenames from input
        self.mummer_filename = mummer_filename
        self.repeats_filename = repeats_filename
        # infer pathlocs and path filenames from mummer filename
        self.pathlocs_filename = mummer_filename.split("mummer")[0] +\
            "pathlocs.txt"
        self.path_filename = mummer_filename.split("mummer")[0] +\
            "paths.txt"
        # infer k (de bruijn graph parameter) and repeat length from filenames
        self.k = mummer_filename.split(".")[1][1:]
        self.min_length = repeats_filename.split(".")[2]
        # get path info from pathlocs and path files
        self.get_path_info()
        # get snp length and termination length in mummer encoding from mummer
        # file
        self.get_mummer_params()
        assert self.termination_length != 1
        # get a long string of the mu mmer-encoded fasta file of paths through
        # SNP graph
        self.get_long_string()
        # get two dictionaries: one from start/end position to path ids and one
        # from start/end position to path names
        self.get_path_locs()
        # self.get_decoded_paths()

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
        for block in self.blocks:
            block.write_to_file(f)
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
        print("Processing mummer file...")
        counter = 0
        for line in lines:
            counter += 1
            if counter % 100000 == 0:
                print("Looking at line {} of mummer file".format(counter))
            start1 = int(line.split()[0]) - 1
            start2 = int(line.split()[1]) - 1
            length = int(line.split()[2])
            self.check_this_match(start1, start2, length)
            og_length = length
            # crop this repeat
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
        for weight in set(weights):
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
        snps/positions."""
        f = open(self.path_filename)
        counter = 0
        self.path_info = dict()
        for line in f:
            if counter == 0:
                # name of path
                name = line.strip()
            elif counter == 1:
                # snps in path
                snps = line.strip().split()
            elif counter == 2:
                # start positions of snps
                positions = [int(x) for x in line.strip().split()]
                self.path_info[name] = (snps, positions)
            counter += 1
            counter = counter % 3

    def get_position(self, name, snps):
        """For a given path and set of snps, return the genetic positions at
        the start and end fo the set of snps, based on the path."""
        return (1, 1)

    def compute_blocks(self, test=False):
        """Put repeats into MPPHB form."""
        self.blocks = []

        for repeat in self.repeats:
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
            # figure out the location on the genome of the first and last snp
            gen_pos1, gen_pos2 = self.get_position(names[0], snps)
            if len(snps) > 0 and len(set(indices)) > 1:
                self.blocks.append(Block(snps, list(set(indices)),
                                         list(set(names))))

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

    def get_decoded_paths(self):
        """Decode paths from long string."""
        long_string_copy = self.mummer_file_string
        paths = []
        while len(long_string_copy) > 0:
            snps = []
            long_string_copy = long_string_copy[1:]
            while long_string_copy[0] != 'X' \
                    and long_string_copy[0] != 'Y':
                # move forward one to account for strt char
                snp = long_string_copy[:self.snp_length]
                long_string_copy = long_string_copy[self.snp_length:]
                snp = snp.replace("C", "0")
                snp = snp.replace("G", "1")
                snp_id = str(int(snp, 2))
                snps.append(snp_id)
                # move forward two to account for SNP call and start char
                long_string_copy = long_string_copy[2:]
            # move forward to next snp
            long_string_copy = long_string_copy[self.termination_length - 1:]
            paths.append(snps)
        self.decoded_paths = paths

    def compute_k(self, snps):
        """For a list of SNPs, compute the number of paths including this
        list."""
        count = 0
        for path in self.decoded_paths:
            if ''.join(snps) in ''.join(path):
                count += 1
        return count

    def compute_selection_coefficients(self):
        """Compute the selection coefficient for each block."""
        # print("Computing k values...")
        for block in self.blocks:
            snps = block.snps
            snps = [x[:-2] for x in snps]
            # print("--Processing block with snps {}".format(snps))
            k = self.compute_k(snps)
            # print("k=", k)
            block.compute_selection_coefficient(k)
