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

def process_repeats(repeats, long_string, locs):
    for repeat in repeats:
        print("Repeat of length {} with {} occurrences".format(repeat[1],
                                    len(repeat[0])))
        starts = repeat[0]
        length = repeat[1]
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
        # figure out how many snps
        start_index = occ.find('A')
        if occ[start_index + 1] == 'A':
            start_index += 1
        overall_length = len(occ) - start_index
        snps = overall_length//SNP_LENGTH
        print("Number of snps: {}".format(snps))
        # figure out which paths
        closest_to_end = 10000000000000
        for start in starts:
            end = start + length - 1
            for key in locs:
                if start >= key[0] and end <= key[1]:
                    name = locs[key]
            print(name)


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

if __name__ == "__main__":
    SNP_LENGTH = 18
    FILENAME =  "yeast10_k100.fa"
    TERMINATION_LENGTH = 13

    f = open(FILENAME, "r")
    # get rid of header line
    f.readline()
    lines = f.readlines()
    lines = [x.strip() for x in lines]
    long_string = ''.join(lines).strip()

    locs = get_path_locs(TERMINATION_LENGTH)

    clean_repeats = get_distinct_repeats(long_string)
    process_repeats(clean_repeats, long_string, locs)
