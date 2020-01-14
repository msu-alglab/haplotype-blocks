def compute_stats(filename):
    """Given a pangenome hap output file, compute avg paths and avg snps"""
    f = open(filename)
    x = []
    y = []
    counter = 0
    for line in f.readlines():
        if counter == 0:
            x.append(len(line.split()))
            counter += 1
        else:
            if line != "No SNPS\n":
                y.append(len(line.split()))
            else:
                print("no SNPS")
                y.append(0)
            counter += 1
        counter = counter % 2
    k = filename.split(".")[1][1:]
    print("For k= {}".format(k))
    total_paths = sum(x)
    avg_paths = total_paths/len(x)
    print("Avg paths:", avg_paths)
    total_snps = sum(y)
    avg_snps = total_snps/len(y)
    print("Avg snps:", avg_snps)


if __name__ == "__main__":
    for filename in [
            "output.k1000.20.txt",
            "output.k100.50.txt",
            "output.k25.20.txt"
    ]:
        compute_stats(filename)
