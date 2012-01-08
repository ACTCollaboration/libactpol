//
// actpol/array.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "state.h"
#include "quaternion.h"

typedef struct
{
    Quaternion focalplane_q;
}
ACTpolFeedhorn;

void
ACTpolFeedhorn_init(ACTpolFeedhorn *feedhorn, double focalplane_x,
    double focalplane_y, double pol_angle);

typedef struct
{
    ACTpolFeedhorn *horn;
    int nhorns;
    double freq_GHz;
}
ACTpolArray;

ACTpolArray *
ACTpolArray_alloc(int nhorns);

void
ACTpolArray_free(ACTpolArray *array);

void
ACTpolArray_init(ACTpolArray *array, double freq_GHz);

typedef struct
{
    const ACTpolArray *array;
    double *ref;
    double mean_ref;
    double *ra, *dec;
    double *cos2gamma, *sin2gamma;
}
ACTpolArrayCoords;

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array);

void
ACTpolArrayCoords_free(ACTpolArrayCoords *coords);

void
ACTpolArrayCoords_init(ACTpolArrayCoords *coords);

// refraction correction only needs to be calculated once per scan,
// or when the weather changes.
void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords,
    const ACTpolScan *scan, const ACTpolWeather *weather);

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state);

int
ACTpolArrayCoords_update_fast(ACTpolArrayCoords *coords, const ACTpolState *state);

#ifdef __cplusplus
}
#endif

