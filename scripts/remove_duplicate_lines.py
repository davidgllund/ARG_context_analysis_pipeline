#!/usr/bin/env python
import sys

file=open(sys.argv[1], 'r')
data = file.readlines()
file.close()

unique=set(data)

output=""
for row in unique:
    output+=str(row)

outfile=open(sys.argv[2], 'w')
outfile.write(output)
outfile.close()
