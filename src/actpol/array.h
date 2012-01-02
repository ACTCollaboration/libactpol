//
// actpol/array.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "state.h"

typedef struct
{
    double freq_GHz;
    double boresight_offset_alt;
    double boresight_offset_az;
    int nhorns;
    double *alt_offset;
    double *az_offset;
}
ACTpolArray;

ACTpolArray *
ACTpolArray_alloc(int ndets);

void
ACTpolArray_free(ACTpolArray *array);

int
ACTpolArray_detector_alt_az(const ACTpolArray *array, int index,
        const ACTpolState *state, double *alt, double *az);

typedef struct
{
    const ACTpolArray *array;
    double *ref;
    double *az, *alt;
    double *ra, *dec;
    double *cos2alpha, *sin2alpha;
}
ACTpolArraySnapshot;

ACTpolArraySnapshot *
ACTpolArraySnapshot_alloc(const ACTpolArray *array);

void
ACTpolArraySnapshot_free(ACTpolArraySnapshot *array);

void
ACTpolArraySnapshot_update_refraction(ACTpolArraySnapshot *snapshot,
        const ACTpolState *state);

int
ACTpolArraySnapshot_update_coords(ACTpolArraySnapshot *snapshot,
        const ACTpolState *state);

#ifdef __cplusplus
}
#endif

