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

    def __init__(self, snps, paths, pathnames, r_method, y_0,
                 length=None, genome_length=None, freqs=None):
        """Inputs:
        snps is a list of snps involved in the block (including 0/1). paths is
        a list of path ids in the block.  pathnames is a list of path names in
        the block. r_method is a string indicating which method should be used
        to compute r.
        Optional parameters:
        for computing the recombination ratio using the length of the block
        relative to the length of the genome, length is the length of this
        block in base pairs, and genome_length is the length of the organism's
        genome in base pairs.
        for computing recombination ratio using frequencies from a data set,
        input a list of four frequencies for the start and end SNPs of the
        block: frequency of start 0 followed by end 0, start 0 followed by end
        1, start 1 followed by end 0, and start 1 followed by end 1."""
        self.snps = snps
        self.paths = paths
        self.pathnames = pathnames
        self.length = length
        self.genome_length = genome_length
        self.compute_recombination_frequency(r_method, freqs)
        self.y_0 = y_0

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
        num_s_to_test = 10000
        self.pop_size = k  # save the k value passed in
        s_to_test = [10*x/num_s_to_test for x in range(num_s_to_test + 1)][1:]
        self.y_t = len(self.paths) / k
        max_likelihood = -1
        best_s = -1
        # print("Computing selection coefficient for block.")
        if self.y_t < 1:
            for s in s_to_test:
                mult1 = 1/s
                to_log = (self.y_t*(1 - self.y_0))/(self.y_0*(1 - self.y_t))
                if to_log > 0:
                    # if not, just do nothing
                    mult2 = math.log(to_log)
                    t = mult1 * mult2
                    summand_1 = -self.r * t
                    to_log = 1-self.y_0*(1-math.e**(s*t))
                    if to_log > 0:
                        # if not, do nothing
                        summand_2 = (self.r/s)*math.log(to_log)
                        inner_log = 1-self.y_0*(1-math.e**(s*t))
                        if inner_log > 0:
                            to_log = t - (1/s) * math.log(inner_log)
                            if to_log > 0:
                                summand_3 = math.log(to_log)
                                likelihood = summand_1 + summand_2 + summand_3
                                if likelihood > max_likelihood:
                                    max_likelihood = likelihood
                                    best_s = s
            self.selection_coefficient = best_s
        else:
            self.selection_coefficient = 0
            # 0 if y_1=1, -1 if logs didn't work out

    def compute_recombination_frequency(self, method=None, freqs=None):
        """Compute the recombination frequency for this block.
        If method is "frequency", need to pass in frequencies of start and end
        of block occurring together."""
        # if method is None, we assume there is no recombination and set r to 0
        if method is None:
            # no recombination
            self.r = 0
        elif method == "distance":
            # compute recombination factor based on number nucleotide between
            # start and end of this block.
            # TODO: choose principled denominator value
            self.r = (1 - math.exp(-self.length/20000)) / 2
        elif method == "frequency":
            # compute recombination fraction based on the data given
            numerator = min(freqs[0] + freqs[3], freqs[1] + freqs[2])
            self.r = numerator / sum(freqs)
            # print(self.r)
            assert self.r >= 0
            assert self.r <= 0.5
