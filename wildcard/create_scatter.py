import numpy as np
import matplotlib.pyplot as plt

files = ["/home/bmumey/wh/chr22-full/output.txt.dist-0.0-500000.txt"]

for file in files:
    my_data = np.genfromtxt(file, delimiter=",")
    sampled = my_data[np.random.choice(my_data.shape[0], 10000, replace=False)]
    x = sampled[:, 0]
    y = sampled[:, 1]
    plt.scatter(x, y, color="red", label="500000")
    plt.xlabel("Number of paths")
    plt.ylabel("Number of SNPs")

    plt.legend()
    plt.savefig("scatterplot.pdf")
