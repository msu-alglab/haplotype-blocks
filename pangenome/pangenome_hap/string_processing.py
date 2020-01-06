import math
import sys

def create_string(filename, output_filename):
    f = open(filename, "r")
    lines = f.readlines()
    num_sequences = len(lines)
    g = open(output_filename, "w")
    counter = 0
    width = math.ceil(math.log2(num_sequences))
    print("width is {}".format(width))
    print("num_sequences is {}".format(num_sequences))
    g.write(">header\n")
    # assume first line is a header, so don't write.
    for index in range(1, num_sequences):
        # look at every line. if header, do one thing, if not, do another.
        print(lines[index][0])
        if lines[index][0] == ">":
            termination_char = format(counter, 'b').zfill(width)
            termination_char = termination_char.replace('0', 'x')
            termination_char = termination_char.replace('1', 'y')
            g.write(termination_char + '\n')
            counter += 1
        else:
            g.write(lines[index])

if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    create_string(input, output)

