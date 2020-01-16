import numpy as np
import matplotlib.pyplot as plt

files = ["/home/bmumey/wh/chr22-full/output.txt.dist-0.0-500000.txt",
         "/home/bmumey/wh/chr22-full/output.txt.dist-0.05-500000.txt",
         "/home/bmumey/wh/chr22-full/output.txt.dist-0.1-500000.txt"]
colors = ["red", "orange", "blue"]

fig, axes = plt.subplots(1, 3)
fig.set_size_inches(10, 4)

ax_index = 0

for file, color in zip(files, colors):
    my_data = np.genfromtxt(file, delimiter=",")
    sampled = my_data[np.random.choice(my_data.shape[0], 10000, replace=False)]
    x = sampled[:, 0]
    y = sampled[:, 1]
    prop = file.split("-")[2]
    ax = axes[ax_index]
    ax.scatter(x, y, color=color, alpha=0.05)
    ax.set_title(prop + " Proportion Wildcard")
    ax_index += 1
    if ax_index == 1:
        ax.set(ylabel="Number of paths")
    ax.set(xlabel="Number of SNPs")

plt.savefig("three_scatterplots.pdf")
