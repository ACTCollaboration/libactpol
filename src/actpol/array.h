//
// actpol/array.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "state.h"
#include "quaternion.h"

#ifdef __cplusplus
extern "C" {
#endif

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
    Quaternion focalplane_q;
}
ACTpolArray;

ACTpolArray *
ACTpolArray_alloc(int nhorns);

void
ACTpolArray_free(ACTpolArray *array);

void
ACTpolArray_init(ACTpolArray *array, double freq_GHz, double focalplane_x, double focalplane_y);

ACTpolFeedhorn *
ACTpolArray_get_feedhorn(ACTpolArray *array, int i);

typedef struct
{
    double a, b;
    double sin2gamma, cos2gamma;
}
ACTpolFeedhornCoords;

enum ACTpolCoordinateSystem {
    ACTPOL_COORDSYS_RA_DEC,
    ACTPOL_COORDSYS_RA_SINDEC,
    ACTPOL_COORDSYS_AZ_ALT,
    ACTPOL_COORDSYS_GALACTIC
};

typedef struct
{
    const ACTpolArray *array;
    double *ref;
    double mean_ref;
    enum ACTpolCoordinateSystem coordsys;
    ACTpolFeedhornCoords *horn;
}
ACTpolArrayCoords;

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array);

void
ACTpolArrayCoords_free(ACTpolArrayCoords *coords);

void
ACTpolArrayCoords_init(ACTpolArrayCoords *coords, enum ACTpolCoordinateSystem coordsys);

// refraction correction only needs to be calculated once per scan,
// or when the weather changes.
void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords,
    const ACTpolScan *scan, const ACTpolWeather *weather);

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state);

ACTpolFeedhornCoords *
ACTpolArrayCoords_get_feedhorn_coords(ACTpolArrayCoords *coords, int i);

#ifdef __cplusplus
}
#endif

