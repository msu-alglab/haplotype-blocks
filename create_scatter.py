import matplotlib.pyplot as plt

f = open("output.txt")
x = []
y = []
counter = 0
for line in f.readlines():
    if counter == 0:
        x.append(len(line.split()))
        counter += 1
    else:
        y.append(len(line.split()))
        counter += 1
    counter = counter % 2

plt.scatter(x,y)
plt.savefig("scatterplot.png")
