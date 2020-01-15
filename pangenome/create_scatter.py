import matplotlib.pyplot as plt

f = open("output.k100.19.txt")
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
            y.append(0)
            print("no SNPS")
        counter += 1

    counter = counter % 2

plt.scatter(x, y, c="red", alpha=0.1, label="k=100")

f = open("output.k500.19.txt")
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

plt.scatter(x, y, c="blue", alpha=0.3, label="k=500")

plt.xlabel("Number of paths")
plt.ylabel("Number of SNPs")

plt.legend()
plt.savefig("scatterplot.pdf")
