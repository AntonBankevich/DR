#!/usr/bin/python3
import sys
aligned = {}
snps = {}
starts = []
finishes = []
#<coord.filtered> <snps>
for line in open (sys.argv[1]):
    arr = line.split()
    starts.append(int(arr[0]))
    finishes.append(int(arr[1]))
starts.sort()
finishes.sort()
print (starts)
print (finishes)
uncovered = 0
total = 0
for line in open (sys.argv[2]):
    arr = line.split()
    coord = int(arr[2])    
    ind = 0
    covered = False
    while ind < len (starts):
        if finishes[ind] < coord:
            ind+=1
            continue
        elif  starts[ind] <= coord:
            covered = True
            break
        else:
            break
    if not covered:
        uncovered +=1
    total +=1
print (f'Total {total} snps, uncovered: {uncovered}')
