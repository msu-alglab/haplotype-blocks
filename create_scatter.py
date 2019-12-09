import matplotlib.pyplot as plt


f = open("output.k100.50.txt")
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
        counter += 1

    counter = counter % 2
indices = [i for (i, val) in enumerate(y) if val != 0]
x = [v for (i, v) in enumerate(x) if i in indices]
y = [v for (i, v) in enumerate(y) if i in indices]

print(indices)
plt.scatter(x,y, c="red", alpha=0.1)

f = open("output.k1000.20.txt")
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
        counter += 1
    counter = counter % 2

indices = [i for (i, val) in enumerate(y) if val != 0]
x = [v for (i, v) in enumerate(x) if i in indices]
y = [v for (i, v) in enumerate(y) if i in indices]
plt.scatter(x,y, c="blue", alpha=0.3)

plt.xlabel("Number of paths")
plt.ylabel("Number of SNPs")
plt.savefig("scatterplot.png")
