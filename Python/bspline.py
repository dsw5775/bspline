def Mn(n, u):
    if n < 2:
        print("Order of B-spline is less than 2\n")
        exit()
    
    if n == 2:
        if u < 0.0 or u > 2.0:
            return 0.0
        
        if u < 1.0:
            return u
        else:
            return 2.0 - u

    return u/(n-1)*Mn(n-1, u) + (n-u)/(n-1)*Mn(n-1, u-1)

def dMn(n, u):
    return Mn(n-1, u) - Mn(n-1, u-1)

def dMn2(n, u):
    return Mn(n-2, u) - 2*Mn(n-2, u-1) + Mn(n-2, u-2)

def bspline_1d(data, min_x1, max_x1, \
    num_x1, x1, order, f_d1 = True, f_d2 = True):

# num_x1:  number of data points
# min_x1:  x-coordinate of the first data point
# max_x1:  x-coordinate of the last data point
# order:   B-spline order
# f_d1:    flag for calculating first-derivative value
# f_d2:    flag for calculating second-derivative value

# value:   returned fitting value
# d1:      returned first-derivative value
# d2:      returned second-derivative value

    numbin_x1 = num_x1 - 1

    x1scaled = (x1 - min_x1) / (max_x1 - min_x1) * numbin_x1

    x1grid0 = int(x1scaled - order/2.0) + 1
    x1grid1 = int(x1scaled + order/2.0)

    if x1grid1 < 0 or x1grid0 > numbin_x1:
        print("out of range (bspline_1d)")
        exit()
    
# If index is nearby the boundary, adjust accordingly.
    if x1grid1 > numbin_x1:
        x1grid1 = numbin_x1
    if x1grid0 < 0:
        x1grid0 = 0

# Mn(order, ...) has domain [0, order] with center (peak) at order/2.0,
# so later on the calculated relative distance to the grid point (x1scaled - i)
# need to be shifted to correspond to the peak.
    x1scaled += order/2.0

    value = 0.0
    d1 = 0.0
    d2 = 0.0    

    for i in range(x1grid0, x1grid1 + 1):
        datapoint = data[i]

        value += datapoint * Mn(order, x1scaled - i)

        if f_d1:
            d1 += datapoint * dMn(order, x1scaled - i) \
                * (numbin_x1 / (max_x1 - min_x1))

        if f_d2:
            d2 += datapoint * dMn2(order, x1scaled - i) \
                * (numbin_x1 / (max_x1 - min_x1)) * (numbin_x1 / (max_x1 - min_x1));

    if f_d2:
        return value, d1, d2
    elif f_d1:
        return value, d1
    else:
        return value

def bspline_2d(data, min_x1, max_x1, num_x1, \
    min_x2, max_x2, num_x2, \
    x1, x2, order, f_d = True):
# f_d:    flag for calculating partial first-derivative value
 
# value:  returned fitting value
# d1:     returned partial first-derivative value along x1 direction
# d2:     returned partial first-derivative value along x2 direction (not to be confused with the one in bspline_1d)

    numbin_x1 = num_x1 - 1
    numbin_x2 = num_x2 - 1

    x1scaled = (x1 - min_x1) / (max_x1 - min_x1) * numbin_x1
    x2scaled = (x2 - min_x2) / (max_x2 - min_x2) * numbin_x2

    x1grid0 = int(x1scaled - order/2.0) + 1
    x1grid1 = int(x1scaled + order/2.0)

    x2grid0 = int(x2scaled - order/2.0) + 1
    x2grid1 = int(x2scaled + order/2.0)

    if x1grid1 < 0 or x1grid0 > numbin_x1 or x2grid1 < 0 or x2grid0 > numbin_x2:
        print("out of range (bspline_2d)")
        exit()
    
    if x1grid1 > numbin_x1:
        x1grid1 = numbin_x1
    if x1grid0 < 0:
        x1grid0 = 0

    if x2grid1 > numbin_x2:
        x2grid1 = numbin_x2
    if x2grid0 < 0:
        x2grid0 = 0

    x1scaled += order/2.0
    x2scaled += order/2.0

    value = 0.0
    d1 = 0.0
    d2 = 0.0    

    for i in range(x1grid0, x1grid1 + 1):
        for j in range(x2grid0, x2grid1 + 1):
            datapoint = data[i][j]

            Mn1 = Mn(order, x1scaled - i)
            Mn2 = Mn(order, x2scaled - j)

            value += datapoint * Mn1 * Mn2

            if f_d:
                dMn1 = dMn(order, x1scaled - i)
                dMn2 = dMn(order, x2scaled - j)

                d1 += datapoint * dMn1 * Mn2 \
                    *(numbin_x1/(max_x1-min_x1))

                d2 += datapoint * Mn1 * dMn2 \
                    *(numbin_x2/(max_x2-min_x2))


    if f_d:
        return value, d1, d2
    else:
        return value