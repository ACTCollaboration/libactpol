/*
 * actpol/array.h : libactpol header file
 *
 * 2011 Mike Nolta <mike@nolta.net>
 */

#pragma once

#include "quaternion.h"

typedef struct
{
    int ndets;
    double boresight_offset_alt;
    double boresight_offset_az;
    double *alt_offset;
    double *az_offset;
}
ACTpolArray;

typedef struct
{
    ACTpolArray *array;
    double *az, *alt;
    double *ra, *dec;
    double *cos2alpha, *sin2alpha;
}
ACTpolArraySnapshot;

ACTpolArray *
ACTpolArray_alloc(int ndets);

void
ACTpolArray_free(ACTpolArray *array);

int
actpol_get_detector_alt_az(const ACTpolArray *array, int index,
        double boresite_alt, double boresite_az, double *alt, double *az);

int
ACTpolArray_compute_snapshot(ACTpolArray *array, Quaternion q, ACTpolArraySnapshot *snapshot);

#ifdef __cplusplus
}
#endif

