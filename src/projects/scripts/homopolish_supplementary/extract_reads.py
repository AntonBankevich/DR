#!/usr/bin/python3
import sys
used_names = set()
for line in open (sys.argv[1]):
    id = line.split()[0]
    used_names.add(id)
cur_pos = 0
to_use = False
id = ""
#grep -P 'ccs 21|ccs -21'
for line in open (sys.argv[2]):
    if line[0] == ">" or line[0] == "@":
        id = line.strip().split()[0][1:]
        if id in used_names:
            print (line[0] + id)        
    else:        
        if id in used_names:
            print (line.strip())
