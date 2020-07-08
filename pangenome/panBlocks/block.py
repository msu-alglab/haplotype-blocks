import math


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
