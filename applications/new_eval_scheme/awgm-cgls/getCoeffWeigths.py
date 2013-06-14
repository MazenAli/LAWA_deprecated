#! /opt/local/bin/python

import sys

sys.stdin.readline()
sys.stdin.readline()

for line in sys.stdin:
    if len(line) == 1:
        continue
        
    line_list = line.split()
    index = line_list[0].split(",")
    type_x = 0 if index[0]=="[scaling" else 1
    type_y = 0 if index[3]=="scaling" else 1
    w = line_list[len(line_list)-1].split("]")
    print type_x, int(index[1]), int(index[2]), type_y, int(index[4]), int(index[5]), float(w[0])    
