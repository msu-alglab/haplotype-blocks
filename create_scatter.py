import matplotlib.pyplot as plt

f = open("output.k1000.20.txt")
x = []
y = []
counter = 0
for line in f.readlines():
    if counter == 0:
        x.append(len(line.split()))
        counter += 1
    else:
        if line != "No SNPS":
            y.append(len(line.split()))
        counter += 1
    counter = counter % 2

plt.scatter(x,y, c="blue")

f = open("output.k100.50.txt")
x = []
y = []
counter = 0
for line in f.readlines():
    if counter == 0:
        x.append(len(line.split()))
        counter += 1
    else:
        if line != "No SNPS":
            y.append(len(line.split()))
        counter += 1
    counter = counter % 2

plt.scatter(x,y, c="red")
plt.savefig("scatterplot.png")
