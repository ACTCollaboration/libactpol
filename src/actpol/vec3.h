
#pragma once

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline double
vec3_norm(const double v[3])
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

static inline void
matrix_times_vec3(double a[3], double m[3][3], const double v[3])
{
    a[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
    a[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
    a[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

static inline void
vec3_cross_product(const double a[3], const double b[3], double axb[3])
{
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
}

#ifdef __cplusplus
}
#endif

