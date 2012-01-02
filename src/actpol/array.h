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
ACTpolArray_center_alt_az(const ACTpolArray *array,
        const ACTpolState *state, double *alt, double *az);

int
ACTpolArray_horn_alt_az(const ACTpolArray *array, int index,
        const ACTpolState *state, double *alt, double *az);

typedef struct
{
    const ACTpolArray *array;
    double *ref;
    double *az, *alt;
    double *ra, *dec;
    double *cos2alpha, *sin2alpha;
}
ACTpolArrayCoords;

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array);

void
ACTpolArrayCoords_free(ACTpolArrayCoords *array);

// refraction correction only needs to be calculated once per scan,
// or when the weather changes.
void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords,
        const ACTpolState *state);

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state);

#ifdef __cplusplus
}
#endif

