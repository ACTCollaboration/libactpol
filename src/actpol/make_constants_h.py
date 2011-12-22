#!/usr/bin/env python

import math

HEADER = """#pragma once
"""

def s(x):
    return "%g" % x

def h(x):
    return x.hex()

def main():

    latitude_deg = -22.9585
    longitude_east_deg = -67.7876
    elevation_m = 5188.0

    d2r = math.pi/180
    latitude = latitude_deg*d2r
    longitude_east = longitude_east_deg*d2r

    defines = [
        ("ACTPOL_LATITUDE_DEG", s(latitude_deg)),
        ("ACTPOL_LONGITUDE_EAST_DEG", s(longitude_east_deg)),
        ("ACTPOL_ELEVATION_METERS", s(elevation_m)+"."),
        ("ACTPOL_LATITUDE", h(latitude)),
        ("ACTPOL_LONGITUDE_EAST", h(longitude_east)),
        ("ACTPOL_SIN_LATITUDE", h(math.sin(latitude))),
        ("ACTPOL_COS_LATITUDE", h(math.cos(latitude))),
        ("ACTPOL_SIN_LONGITUDE_EAST", h(math.sin(longitude_east))),
        ("ACTPOL_COS_LONGITUDE_EAST", h(math.cos(longitude_east))),
    ]

    f = open("constants.h", "w")
    f.write(HEADER)
    for define in defines:
        #print "%32s  %s" % define
        f.write("#define %s %s\n" % define)
    f.close()

if __name__ == '__main__':
    main()

