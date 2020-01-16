import matplotlib.pyplot as plt

files = ["/home/bmumey/wh/chr22-full/output.txt.dist-0.0-500000.txt"]
colors = ["red"]

for file, color in zip(files, colors):
    f = open(file)
    lines = f.readlines()
    x = [x.split(",")[0] for x in lines]
    y = [x.split(",")[1] for x in lines]
    plt.scatter(x, y, c=color, alpha=0.1, label="k=100")

plt.xlabel("Number of paths")
plt.ylabel("Number of SNPs")

plt.legend()
plt.savefig("scatterplot.pdf")
