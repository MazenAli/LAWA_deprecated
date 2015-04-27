#! /opt/local/bin/python

import sys
count = 0
for line in sys.stdin:
    entries = line.split()
    print " ".join(entries)
    count +=1
    if count==int(sys.argv[1]):
        count = 0
        print
