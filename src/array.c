
#include <assert.h>
#include <stdlib.h>

#include "actpol/array.h"
#include "actpol/astro.h"
#include "actpol/state.h"
#include "actpol/vec3.h"

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

int
ACTpolArray_center_alt_az(const ACTpolArray *array,
        const ACTpolState *state, double *alt, double *az)
{
    *alt = state->boresight_alt
         + array->boresight_offset_alt;
    *az = state->boresight_az
        + array->boresight_offset_az;

    return 0;
}

int
ACTpolArray_horn_alt_az(const ACTpolArray *array, int index,
        const ACTpolState *state, double *alt, double *az)
{
    assert(index >= 0 && index < array->nhorns);

    *alt = state->boresight_alt
         + array->boresight_offset_alt
         + array->horn[index].alt_offset;
    *az = state->boresight_az
        + array->boresight_offset_az
        + array->horn[index].az_offset;

    return 0;
}

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array)
{
    ACTpolArrayCoords *coords = (ACTpolArrayCoords *)malloc(sizeof(ACTpolArrayCoords));
    coords->array = array;
    coords->ref = (double *)malloc(sizeof(double) * array->nhorns);
    coords->ra = (double *)malloc(sizeof(double) * array->nhorns);
    coords->dec = (double *)malloc(sizeof(double) * array->nhorns);
    coords->sin2alpha = (double *)malloc(sizeof(double) * array->nhorns);
    coords->cos2alpha = (double *)malloc(sizeof(double) * array->nhorns);
    return coords;
}

void
ACTpolArrayCoords_free(ACTpolArrayCoords *coords)
{
    free(coords->cos2alpha);
    free(coords->sin2alpha);
    free(coords->dec);
    free(coords->ra);
    free(coords->ref);
    free(coords);
}

void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords, const ACTpolState *state)
{
    const ACTpolArray *array = coords->array;
    for (int i = 0; i != coords->array->nhorns; ++i)
    {
        double alt, az;
        ACTpolArray_horn_alt_az(coords->array, i, state, &alt, &az);
        coords->ref[i] = actpol_refraction(&state->weather, array->freq_GHz, alt);
    }
}

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state)
{
    const ACTpolArray *array = coords->array;
    Quaternion focalplane_to_GCRS;
    Quaternion_mul(focalplane_to_GCRS, state->NWU_to_GCRS_q, state->focalplane_to_NWU_q);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_GCRS, array->horn[i].focalplane_q);

        double mat[3][3];
        Quaternion_to_matrix(q, mat);
        const double r[3] = {mat[0][2], mat[1][2], mat[2][2]};

        actpol_vec2ang(r, coords->ra+i, coords->dec+i);
    }

    return 0;
}

