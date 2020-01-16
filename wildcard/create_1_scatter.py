import numpy as np
import matplotlib.pyplot as plt

files = ["/home/bmumey/wh/chr22-full/output.txt.dist-0.1-500000.txt",
         # "/home/bmumey/wh/chr22-full/output.txt.dist-0.05-500000.txt",
         "/home/bmumey/wh/chr22-full/output.txt.dist-0.0-500000.txt"]
colors = ["blue",
          # "orange",
          "red"]


for file, color in zip(files, colors):
    my_data = np.genfromtxt(file, delimiter=",")
    sampled = my_data[np.random.choice(my_data.shape[0], 10000, replace=False)]
    x = sampled[:, 0]
    y = sampled[:, 1]
    prop = file.split("-")[2]
    plt.scatter(x, y, color=color, alpha=0.05, label=prop, s=10)

plt.xlabel("Number of paths")
plt.ylabel("Number of SNPs")
plt.legend()
plt.savefig("1_scatterplot.pdf")
