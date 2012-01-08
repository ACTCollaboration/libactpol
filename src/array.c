
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "actpol/array.h"
#include "actpol/astro.h"
#include "actpol/state.h"
#include "actpol/vec3.h"

#include "debug.h"

void
ACTpolFeedhorn_init(ACTpolFeedhorn *feedhorn, double focalplane_x,
    double focalplane_y, double pol_angle)
{
    Quaternion_r3(feedhorn->focalplane_q, -pol_angle);
    Quaternion_r2_mul(focalplane_x, feedhorn->focalplane_q);
    Quaternion_r1_mul(-focalplane_y, feedhorn->focalplane_q);
}

ACTpolArray *
ACTpolArray_alloc(int nhorns)
{
    ACTpolArray *array;
    array = (ACTpolArray *)malloc(sizeof(ACTpolArray));
    assert(array);
    array->nhorns = nhorns;
    array->horn = (ACTpolFeedhorn *)malloc(nhorns*sizeof(ACTpolFeedhorn));
    assert(array->horn);
    return array;
}

void
ACTpolArray_free(ACTpolArray *array)
{
    free(array->horn);
    free(array);
}

void
ACTpolArray_init(ACTpolArray *array, double freq_GHz)
{
    array->freq_GHz = freq_GHz;
}

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array)
{
    ACTpolArrayCoords *coords = (ACTpolArrayCoords *)malloc(sizeof(ACTpolArrayCoords));
    coords->array = array;
    coords->ref = (double *)malloc(sizeof(double) * array->nhorns);
    coords->ra = (double *)malloc(sizeof(double) * array->nhorns);
    coords->dec = (double *)malloc(sizeof(double) * array->nhorns);
    coords->sin2gamma = (double *)malloc(sizeof(double) * array->nhorns);
    coords->cos2gamma = (double *)malloc(sizeof(double) * array->nhorns);
    return coords;
}

void
ACTpolArrayCoords_free(ACTpolArrayCoords *coords)
{
    free(coords->cos2gamma);
    free(coords->sin2gamma);
    free(coords->dec);
    free(coords->ra);
    free(coords->ref);
    free(coords);
}

void
ACTpolArrayCoords_init(ACTpolArrayCoords *coords)
{
    coords->mean_ref = arcsec2rad(30.);
}

void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords,
    const ACTpolScan *scan, const ACTpolWeather *weather)
{
    const ACTpolArray *array = coords->array;

    Quaternion focalplane_to_NWU_q;
    Quaternion_identity(focalplane_to_NWU_q);
    actpol_rotate_focalplane_to_NWU(scan->mean_alt, scan->mean_az, focalplane_to_NWU_q);

    coords->mean_ref = 0.;
    for (int i = 0; i != coords->array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_NWU_q, array->horn[i].focalplane_q);

        double mat[3][3];
        Quaternion_to_matrix(q, mat);
        double alt = asin(mat[2][2]);

        coords->ref[i] = actpol_refraction(weather, array->freq_GHz, alt);
        assert(coords->ref[i] > 0. && coords->ref[i] < 1e-3);
        coords->mean_ref += coords->ref[i];
    }
    coords->mean_ref /= coords->array->nhorns;
    DEBUG("mean_ref = %g\"\n", rad2arcsec(coords->mean_ref));
}

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state)
{
    const ACTpolArray *array = coords->array;

    Quaternion focalplane_to_NWU_q;
    Quaternion_identity(focalplane_to_NWU_q);
    actpol_rotate_focalplane_to_NWU(
        state->boresight_alt - coords->mean_ref,
        state->boresight_az,
        focalplane_to_NWU_q);

    Quaternion focalplane_to_GCRS;
    Quaternion_mul(focalplane_to_GCRS, state->NWU_to_GCRS_q, focalplane_to_NWU_q);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_GCRS, array->horn[i].focalplane_q);

        double mat[3][3];
        Quaternion_conj(q);
        Quaternion_to_matrix(q, mat);
        double *p1 = mat[0];
        double *p2 = mat[1];
        double *r = mat[2];

        // XXX:speed this up
        actpol_vec2ang(r, coords->ra+i, coords->dec+i);

        // w = r x z
        double w[3], z[3] = {0, 0, 1};
        vec3_cross_product(w, r, z);
        vec3_unit(w);

        // n = w x r
        double n[3];
        vec3_cross_product(n, w, r);

        double sin_g = vec3_dot_product(p1, w);
        double cos_g = vec3_dot_product(p1, n);
        coords->sin2gamma[i] = 2*sin_g*cos_g;
        coords->cos2gamma[i] = 2*cos_g*cos_g - 1;
    }

    return 0;
}

