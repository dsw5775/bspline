#!/usr/bin/python

from bspline import *
import numpy as np
import time
import sys

def main():
# playground for testing bspline_1d to bspline_2d with order 6
# you can change order with command line parameter, e.g. ./testbspline 4 

    order = 6

    if(len(sys.argv) > 1):
        order = int(sys.argv[1])
        
    num_x1 = 100 + 1

    min_x1 = 0.0
    max_x1 = 1.0
    width_x1 = (max_x1-min_x1) / (num_x1-1)
    x1 = 0.5

    data1 = np.empty((num_x1))

    for i in range(0, num_x1):
        data1[i] = (width_x1*i) * (width_x1*i)

    start = time.perf_counter()
    value, d1, d2 = bspline_1d(data1, min_x1, max_x1, num_x1, x1, order)
    end = time.perf_counter()

    print("bspline_1d test:")
    print("Calculated:")
    print(f'{value:20.6f} {d1:20.6f} {d2:20.6f}')
    print("Expected:")
    print(f'{x1*x1:20.6f} {2*x1:20.6f} {2.0:20.6f}')
    print(f'Time: {end - start:10.6f} sec\n')

# test bspline_2d for data generated from function x1^2*x2
    num_x1 = 100; 
    min_x1 = 0.005
    max_x1 = 0.995
    width_x1 = (max_x1-min_x1) / (num_x1-1)
    x1 = 0.3
    num_x2 = 24+1
    min_x2 = -10000.0
    max_x2 = 10000.0;
    width_x2 = (max_x2-min_x2) / (num_x2-1)
    x2 = 4000.0;

    data2 = np.empty((num_x1, num_x2))

    for i in range(0, num_x1):
        for j in range(0, num_x2):
            data2[i][j] = (min_x1+width_x1*i) * (min_x1+width_x1*i) \
                * (min_x2+width_x2*j)

    start = time.perf_counter()
    value, d1, d2 = bspline_2d(data2, min_x1, max_x1, num_x1, \
        min_x2, max_x2, num_x2, \
        x1, x2, order)
    end = time.perf_counter()

    print("bspline_2d test:")
    print("Calculated:")
    print(f'{value:20.6f} {d1:20.6f} {d2:20.6f}')
    print("Expected:")
    print(f'{x1*x1*x2:20.6f} {2*x1*x2:20.6f} {x1*x1:20.6f}')
    print(f'Time: {end - start:10.6f} sec\n')

if __name__ == "__main__":
    main()
