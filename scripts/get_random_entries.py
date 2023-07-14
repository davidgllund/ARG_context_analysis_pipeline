#!/usr/bin/env python
from random import sample
import sys

file1=open(sys.argv[1], 'r')
head=file1.readlines()
file1.close()

samples=sample(head,int(sys.argv[2]))

stringToExport=""

for row in samples:
    stringToExport += str(row)

fp=open(sys.argv[3], 'w')
fp.write(stringToExport)
fp.close()
