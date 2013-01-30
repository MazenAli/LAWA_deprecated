#! /opt/local/bin/python

import sys

for line in sys.stdin:
    if len(line) == 1:
        continue
        
    line_list = line.split(",")
    type_x = 0 if line_list[0]=="scaling" else 1
    type_y = 0 if line_list[3]=="scaling" else 1
    print type_x, int(line_list[1]), int(line_list[2]), type_y, int(line_list[4]), int(line_list[5])   
