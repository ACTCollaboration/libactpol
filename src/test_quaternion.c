
#include <assert.h>
#include <stdio.h>

#include "actpol/actpol.h"
#include "sofa.h"

void
print_mat(double mat[3][3])
{
    for (int i = 0; i < 3; i++)
        printf("%+.15f %+.15f %+.15f\n", mat[i][0], mat[i][1], mat[i][2]);
}

void
print_quat(Quaternion q)
{
    printf("q: %.15e %.15e %.15e %.15e\n", q[0], q[1], q[2], q[3]);
}

int
main(int argc, char *argv[])
{
    const double dut1 = -0.072073685;
    const double as2r = 4.848136811095359935899141e-6;
    const double xp = 0.0349282 * as2r;
    const double yp = 0.4833163 * as2r;

    double djmjd0, date, time, utc, dat, tai, tt, tut, ut1;
    double X, Y, s, theta;

    iauCal2jd(2007, 4, 5, &djmjd0, &date);
    time = 12./24.;
    utc = date + time;
    iauDat(2007, 4, 5, time, &dat);
    tai = utc + dat/86400.;
    tt = tai + 32.184/86400.;
    tut = time + dut1/86400.;
    ut1 = date + tut;

    iauXy06(djmjd0, tt, &X, &Y);
    s = iauS06(djmjd0, tt, X, Y);
    X += 0.175e-3*as2r;
    Y += -0.2259e-3*as2r;
    printf("X, Y = %e, %e\n", X, Y);
    printf("s = %e\n", s/as2r);

    double rc2i[3][3];
    iauC2ixys(X, Y, s, rc2i);
    print_mat(rc2i);

    double Z = sqrt(1. - X*X - Y*Y);
    double d = acos(Z);
    double E = atan2(Y, X);
    theta = iauEra00(djmjd0+date, tut);
    printf("theta = %.15f\n", theta*180./M_PI);
    printf("Z,E,d = %e %e %e\n", Z, E, d);
    Quaternion q;
    Quaternion_r3(q,-E);
    Quaternion_r2_mul(-d,q);
    Quaternion_r3_mul(E+s,q);
    //Quaternion_r2_mul(xp, q);
    //Quaternion_r1_mul(yp, q);

    double mat[3][3];
    Quaternion_to_matrix(q, mat);
    print_mat(mat);
    return 0;
}

