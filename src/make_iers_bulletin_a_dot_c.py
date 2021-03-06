#!/usr/bin/env python

header = """
#include <assert.h>
#include <math.h>

#include "actpol/iers_bulletin_a.h"

typedef struct { float x; float y; float dut1; } xyz;

"""

footer = """
int
actpol_get_iers_bulletin_a( double mjd, double *dut1, double *x, double *y )
{
    if ( (mjd_min < mjd) && (mjd < mjd_max) )
    {
        double mjd_floor;
        double r = modf( mjd, &mjd_floor );
        int k = (int) mjd_floor - mjd_min;
        xyz a = bulletinA[k];
        xyz b = bulletinA[k+1];

        // Detect leap seconds
        double leap = b.dut1 - a.dut1;
        if (leap > 0.5)
            leap = 1.;
        else if (leap < -0.5)
            leap = -1.;
        else
            leap = 0;

        *dut1 = (1-r)*a.dut1 + r*(b.dut1-leap);
        *x = (1-r)*a.x + r*b.x;
        *y = (1-r)*a.y + r*b.y;
    }
    else
    {
        assert( 1 == 0 );
        return -1;
    }

    return 0;
}

"""

if __name__ == '__main__':

    f = open( "iers_bulletin_a.c", "w" )
    f.write( header )

    mjds = []
    f.write( "static const xyz bulletinA[] = {\n" )

    for line in open('finals2000A.truncated'):
        date = line[:6]
        mjd = int(float(line[7:15]))

        if line[23] != ' ':
            # bulletin A
            x = float(line[18:27])
            y = float(line[37:46])
            dut1 = float(line[58:68])
        else:
            break

        mjds.append(mjd)
        f.write("  {% ff, % ff, % .7ff},  // %s\n" % (x,y,dut1,date))

    f.write( "};\n\n" )
    f.write( "static const int mjd_min = %d;\n" % min(mjds) )
    f.write( "static const int mjd_max = %d;\n" % max(mjds) )

    f.write( footer )
    f.close()

