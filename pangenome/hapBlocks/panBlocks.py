import networkx as nx
import math
import matplotlib.pyplot as plt
import numpy as np


class Block:
    """A Block object contains the SNP nodes (ID and 0/1) and paths (by both ID
    and name) for a maximal pangenome haplotype block. Based on this
    information, it can compute the selection coefficient for each block.

    Attributes:
    * snps: the snps involved in this block
    * paths: the path ids of this block
    * pathnames: the path names of this block
    * length: length along the genome of this block
    * genome_length: total genome length for this species
    """

    def __init__(self, snps, paths, pathnames, length, genome_length):
        self.snps = snps
        self.paths = paths
        self.pathnames = pathnames
        self.length = length
        self.genome_length = genome_length
        self.compute_recombination_frequency(method="distance")
        self.y_0 = 0.00005  # as in Cunha et al.

    def __str__(self):
        """Return a string representation of |K|
        and number snps in this block"""
        return "|K|={}, num snps={}".format(
            len(self.paths),
            len(self.snps)
        )

    def write_to_file(self, f):
        """Write this block's info to passed file."""
        f.write(" ".join(self.paths) + "\n")
        f.write(" ".join(self.snps) + "\n")

    def compute_selection_coefficient(self, k):
        """Compute the selection coefficient for this block.
        Use Equation 5 from Cunha et al."""
        num_s_to_test = 10
        s_to_test = [x/num_s_to_test for x in range(num_s_to_test + 1)][1:]
        y_t = len(self.paths) / k
        max_likelihood = -1
        best_s = -1
        # print("Computing selection coefficient for block.")
        if y_t < 1:
            for s in s_to_test:
                t = (1/s) * math.log((y_t*(1 - self.y_0))/(self.y_0*(1 - y_t)))
                summand_1 = -self.r * t
                summand_2 = (self.r/s)*math.log(1-self.y_0*(1-math.e**(s*t)))
                summand_3 = math.log(t-(1/2)*math.log(
                    -self.y_0*(1-math.e**(s*t))))
                likelihood = summand_1 + summand_2 + summand_3
                if likelihood > max_likelihood:
                    max_likelihood = likelihood
                    best_s = s
            self.selection_coefficient = best_s
            # print("s=", self.selection_coefficient)
        else:
            # print("y_t is 1, so skipping and setting s to 0")
            self.selection_coefficient = 0

    def compute_recombination_frequency(self, method=None):
        """Compute the recombination frequency for this block."""
        # if method is None, we assume there is no recombination and set r to 0
        if method is None:
            # no recombination
            self.r = 0
        elif method == "distance":
            # compute recombination factor based on number nucleotide between
            # start and end of this block.
            # TODO: choose principled denominator value
            self.r = (1 - math.exp(-self.length/20000)) / 2


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
        # TODO: filenames are not consistent. this is broken.
        self.k = mummer_filename.split(".")[1][1:]
        self.min_length = repeats_filename.split(".")[2]
        # get path info from pathlocs and path files
        # also get genome length
        self.get_path_info()
        # set all_paths from path_info
        self.all_paths = [x[0] for x in self.path_info.values()]
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
        snps/positions. Also, create a list of all unique snps."""
        f = open(self.path_filename)
        counter = 0
        overall_max_position = 0
        self.path_info = dict()
        self.snps = dict()
        for line in f:
            if counter == 0:
                # name of path
                name = line.strip()
            elif counter == 1:
                # snps in path
                snps = line.strip().split()
                for snp in snps:
                    self.snps[snp] = 0
            elif counter == 2:
                # start positions of snps
                positions = [int(x) for x in line.strip().split()]
                max_pos = positions[-1]
                if max_pos > overall_max_position:
                    overall_max_position = max_pos
                self.path_info[name] = (snps, positions)
            counter += 1
            counter = counter % 3
        self.genome_length = overall_max_position

    def get_position(self, name, snps):
        """For a given path and set of snps, return the genetic positions at
        the start and end of the set of snps, based on the path."""
        # name is the name of one of the paths with this block. snps is the
        # snps in the block.
        # need to look up the positions from path_info dict
        info = self.path_info[name]
        path_snps = info[0]
        index = path_snps.index(snps[0])
        start_position = info[1][index]
        end_position = info[1][index + len(snps) - 1]
        return (start_position, end_position)

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
            distance = gen_pos2 - gen_pos1
            if len(snps) > 0 and len(set(indices)) > 1:
                self.blocks.append(Block(snps,
                                         list(set(indices)),
                                         list(set(names)),
                                         distance,
                                         self.genome_length
                                         )
                                   )

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
        for block in self.blocks:
            snps = block.snps
            snps = [x[:-2] for x in snps]
            # print("--Processing block with snps {}".format(snps))
            k = self.compute_k(snps)
            block.compute_selection_coefficient(k)

    def plot_snps_vs_paths_scatter(self, filename=None):
        """Create a scatterplot of all blocks in PanBlocks object,
        with number of paths on the x axis and number of snps on the y axis"""
        # if no filename provided, write to wherever the mummer file is
        if filename is None:
            filename = self.mummer_filename + "scatterplot.pdf"
        # x values are number of paths in the block, y are number of snps
        x = np.empty(self.get_num_blocks())
        y = np.empty(self.get_num_blocks())
        for i, block in enumerate(self.get_blocks()):
            x[i] = len(block.snps)
            y[i] = len(block.paths)
        plt.scatter(x, y, c="red", alpha=0.5, label="k={}".format(19))
        plt.xlabel("Number of paths")
        plt.ylabel("Number of SNPs")
        plt.legend()
        # set title to mummer filename. TODO: better title
        plt.title(self.mummer_filename)
        print("Saving scatterplot as {}".format(filename))
        plt.savefig(filename)

    def set_selection_coefficients(self):
        """For each node, figure out what the max selection coefficient is."""
        # TODO: make this more efficient
        for snp in self.snps.keys():
            # get the max selection coeff for this node
            # print(f"Finding max selection coeff for snp {snp}")
            max_s = 0
            for block in self.blocks:
                if snp in block.snps:
                    if block.selection_coefficient > max_s:
                        max_s = block.selection_coefficient
            self.snps[snp] = max_s

    def print_selection_coefficients(self):
        """Print out all selection coeffs."""
        print("Calling print_selection_coefficients")
        print(self.snps)

    def generate_hists(self):
        """Make a histogram of selection coeffs by snp and by block"""
        plt.subplot(211)
        plt.hist(self.snps.values())
        plt.subplot(212)
        # get all selection coeffs from blocks
        coeffs = []
        for block in self.blocks:
            coeffs.append(block.selection_coefficient)
        plt.hist(coeffs)
        plt.savefig("snp_selection_coeffs.png")

    def decorate_snp_graph(self, dotfile):
        f = open(dotfile)
        out = open(dotfile.split(".")[0] + "_new.dot", "w")
        print("Editing dotfile")
        print(self.snps)
        # add in a color based on selection coeff
        counter = 0  # assume that the snps are in order in the dotfile
        one = False
        for line in f:
            if "style=dotted" in line:
                parts = line.split("dotted")
                side = line.split("label=")[1][1]
                snp_label = str(counter) + ":" + side
                selection_coeff = self.snps[snp_label]
                if snp_label == "278:1":
                    print("SNP {} has selection coeff {}".format(
                        snp_label,
                        selection_coeff))
                red_val = f"{int(255 - selection_coeff * 255):0>2x}"
                color = '"#ff' + red_val + red_val + '"'
                new_line = parts[0] + "filled fillcolor=" +\
                    color + " color=black" + parts[1]
                out.write(new_line)
                print(f"line={line}, snp_label={snp_label}, counter={counter}")
                print(color)
                if one:
                    counter += 1
                one = not one
            else:
                out.write(line)
        f.close()
        out.close()
