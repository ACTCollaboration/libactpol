#!/usr/bin/env python

import math

HEADER = """#pragma once
"""

def d(x):
    return repr(x)

def main():

    latitude_deg = -22.9585
    longitude_east_deg = -67.7876
    elevation_m = 5188.0
    speed_of_light = 299792458.

    d2r = math.pi/180
    latitude = latitude_deg*d2r
    longitude_east = longitude_east_deg*d2r

    defines = [
        ("ACTPOL_LATITUDE_DEG", d(latitude_deg)),
        ("ACTPOL_LONGITUDE_EAST_DEG", d(longitude_east_deg)),
        ("ACTPOL_ELEVATION_METERS", d(elevation_m)),
        ("SPEED_OF_LIGHT_M_PER_S", d(speed_of_light)),
        ("ACTPOL_LATITUDE", d(latitude)),
        ("ACTPOL_LONGITUDE_EAST", d(longitude_east)),
        ("ACTPOL_SIN_LATITUDE", d(math.sin(latitude))),
        ("ACTPOL_COS_LATITUDE", d(math.cos(latitude))),
        ("ACTPOL_SIN_LONGITUDE_EAST", d(math.sin(longitude_east))),
        ("ACTPOL_COS_LONGITUDE_EAST", d(math.cos(longitude_east))),
    ]

    f = open("constants.h", "w")
    f.write(HEADER)
    for define in defines:
        #print "%32s  %s" % define
        f.write("#define %s %s\n" % define)
    f.close()

if __name__ == '__main__':
    main()

