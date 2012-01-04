
#include <assert.h>
#include <stdlib.h>

#include "actpol/array.h"
#include "actpol/astro.h"
#include "actpol/state.h"
#include "actpol/vec3.h"

ACTpolArray *
ACTpolArray_alloc(int nhorns)
{
    ACTpolArray *array;
    array = (ACTpolArray *)malloc(sizeof(ACTpolArray));
    assert(array);
    array->nhorns = nhorns;
    array->horn = (ACTpolFeedhorn *)malloc(nhorns*sizeof(ACTpolFeedhorn));
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
    coords->az = (double *)malloc(sizeof(double) * array->nhorns);
    coords->alt = (double *)malloc(sizeof(double) * array->nhorns);
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
    free(coords->alt);
    free(coords->az);
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
    double mat[3][3];
    Quaternion_to_matrix(state->q, mat);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        double alt, az;
        ACTpolArray_horn_alt_az(array, i, state, &alt, &az);
        alt -= coords->ref[i]; // correct for refraction

        coords->alt[i] = alt;
        coords->az[i] = az;

        // rotate horizon -> celestial
        double r_h[3], r_c[3];

        actpol_ang2vec(-az, alt, r_h);
        /*
        double cos_d = cos(alt);
        r_h[0] = cos_d*cos(-az);
        r_h[1] = cos_d*sin(-az);
        r_h[2] = sin(alt);
        */

        matrix_times_vec3(r_c, mat, r_h);

        actpol_vec2ang(r_c, coords->ra+i, coords->dec+i);
        /*
        coords->ra[i] = atan2(r_c[1], r_c[0]);
        coords->dec[i] = atan2(r_c[2], hypot(r_c[0],r_c[1]));
        */
    }

    return 0;
}

