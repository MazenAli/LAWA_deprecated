#! /opt/local/bin/python

import sys

sys.stdin.readline()

for line in sys.stdin:
    if len(line) == 1:
        continue
        
    line_list = line.split(",")
    dummy_list = line_list[0].split("[")
    line_list[0] = dummy_list[1]
    dummy_list = line_list[5].split("]")
    line_list[5] = dummy_list[0]
    type_x = 0 if line_list[0]=="scaling" else 1
    type_y = 0 if line_list[3]=="scaling" else 1
    print type_x, int(line_list[1]), int(line_list[2]), type_y, int(line_list[4]), int(line_list[5])   
