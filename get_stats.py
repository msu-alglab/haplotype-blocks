import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': '14'})

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
            print("no SNPS")
        counter += 1

    counter = counter % 2

#idx = [i for (i,v) in enumerate(y) if v >= 2]
#x = [v for (i,v) in enumerate(x) if i in idx]
#y = [v for (i,v) in enumerate(y) if i in idx]

print("For k=100")
total_paths = sum(x)
avg_paths = total_paths/len(x)
print("Avg paths:", avg_paths)
total_snps = sum(y)
avg_snps = total_snps/len(y)
print("Avg snps:", avg_snps)

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
            print("no SNPS")
            y.append(0)
        counter += 1
    counter = counter % 2
#idx = [i for (i,v) in enumerate(y) if v >= 2]
#x = [v for (i,v) in enumerate(x) if i in idx]
#y = [v for (i,v) in enumerate(y) if i in idx]
print()
print("For k=1000")
total_paths = sum(x)
avg_paths = total_paths/len(x)
print("Avg paths:", avg_paths)
total_snps = sum(y)
avg_snps = total_snps/len(y)
print("Avg snps:", avg_snps)
