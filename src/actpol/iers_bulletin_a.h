
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int
actpol_get_iers_bulletin_a( double mjd, double *dut1, double *x, double *y );

int
actpol_set_iers_bulletin_a( int mjd_min_, int mjd_max_, double *dut1, double *x, double *y);

#ifdef __cplusplus
}
#endif

