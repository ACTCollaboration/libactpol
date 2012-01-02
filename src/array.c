
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
    array->nhorns = nhorns;
    array->alt_offset = (double *)malloc(nhorns*sizeof(double));
    array->az_offset = (double *)malloc(nhorns*sizeof(double));
    return array;
}

void
ACTpolArray_free(ACTpolArray *array)
{
    free(array->az_offset);
    free(array->alt_offset);
    free(array);
}

int
ACTpolArray_center_altaz(const ACTpolArray *array,
        const ACTpolState *state, double *alt, double *az)
{
    *alt = state->boresight_alt
         + array->boresight_offset_alt;
    *az = state->boresight_az
        + array->boresight_offset_az;

    return 0;
}

int
ACTpolArray_detector_alt_az(const ACTpolArray *array, int index,
        const ACTpolState *state, double *alt, double *az)
{
    assert(index >= 0 && index < array->nhorns);

    *alt = state->boresight_alt
         + array->boresight_offset_alt
         + array->alt_offset[index];
    *az = state->boresight_az
        + array->boresight_offset_az
        + array->az_offset[index];

    return 0;
}

ACTpolArraySnapshot *
ACTpolArraySnapshot_alloc(const ACTpolArray *array)
{
    ACTpolArraySnapshot *snapshot = (ACTpolArraySnapshot *)malloc(sizeof(ACTpolArraySnapshot));
    snapshot->array = array;
    snapshot->ref = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->az = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->alt = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->ra = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->dec = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->sin2alpha = (double *)malloc(sizeof(double) * array->nhorns);
    snapshot->cos2alpha = (double *)malloc(sizeof(double) * array->nhorns);
    return snapshot;
}

void
ACTpolArraySnapshot_free(ACTpolArraySnapshot *snapshot)
{
    free(snapshot->cos2alpha);
    free(snapshot->sin2alpha);
    free(snapshot->dec);
    free(snapshot->ra);
    free(snapshot->alt);
    free(snapshot->az);
    free(snapshot->ref);
    free(snapshot);
}

void
ACTpolArraySnapshot_update_refraction(ACTpolArraySnapshot *snapshot, const ACTpolState *state)
{
    const ACTpolArray *array = snapshot->array;
    for (int i = 0; i != snapshot->array->nhorns; ++i)
    {
        double alt, az;
        ACTpolArray_detector_alt_az(snapshot->array, i, state, &alt, &az);
        snapshot->ref[i] = actpol_refraction(&state->weather, array->freq_GHz, alt);
    }
}

int
ACTpolArraySnapshot_update_coords(ACTpolArraySnapshot *snapshot, const ACTpolState *state, const Quaternion q)
{
    const ACTpolArray *array = snapshot->array;
    double mat[3][3];
    Quaternion_to_matrix(q, mat);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        double alt, az;
        ACTpolArray_detector_alt_az(array, i, state, &alt, &az);
        alt -= snapshot->ref[i]; // correct for refraction

        snapshot->alt[i] = alt;
        snapshot->az[i] = az;

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

        actpol_vec2ang(r_c, snapshot->ra+i, snapshot->dec+i);
        /*
        snapshot->ra[i] = atan2(r_c[1], r_c[0]);
        snapshot->dec[i] = atan2(r_c[2], hypot(r_c[0],r_c[1]));
        */
    }

    return 0;
}

#if 0
// atan x = (x/(1+x^2)) (1 + (2/3)z + (8/15)z^2 + ...)
// where z = x^2/(1+x^2)
//
double
rel_atan2(double y0, double x0, double atan0, double y, double x)
{
    double z = y/x;
    //double z0 = y0/x0;
    double atan1 = x0*(z*x0 - y0)/(x0*x0 + y0*y0);
    return atan0 + atan1;
}

int
ACTpolArraySnapshot_update_coords_fast(ACTpolArraySnapshot *snapshot, const ACTpolState *state)
{
    const ACTpolArray *array = snapshot->array;
    double mat[3][3];
    double alt, az, center_ra, center_dec;
    double r_h[3], r_c[3];
    double center_ra_y, center_ra_x;
    double center_dec_y, center_dec_x;

    Quaternion_to_matrix(state->q, mat);

    ACTpolArray_center_altaz(array, state, &alt, &az);
    actpol_ang2vec(-az, alt, r_h);
    matrix_times_vec3(r_c, mat, r_h);
    actpol_vec2ang(r_c, &center_ra, &center_dec);
    center_ra_y = r_h[1];
    center_ra_x = r_h[0];
    center_dec_y = r_h[2];
    center_dec_x = hypot(r_h[0],r_h[1]);

    for (int i = 0; i != array->nhorns; ++i)
    {
        ACTpolArray_detector_alt_az(array, i, state, &alt, &az);
        alt -= snapshot->ref[i]; // correct for refraction
        snapshot->alt[i] = alt;
        snapshot->az[i] = az;

        // rotate horizon -> celestial
        actpol_ang2vec(-az, alt, r_h);
        matrix_times_vec3(r_c, mat, r_h);
        //actpol_vec2ang(r_c, snapshot->ra+i, snapshot->dec+i);
        snapshot->ra[i] = rel_atan2(center_ra_y, center_ra_x, center_ra, r_h[1], r_h[0]);
        snapshot->dec[i] = rel_atan2(center_dec_y, center_dec_y, center_dec, r_h[2], hypot(r_h[0],r_h[1]));
    }

    return 0;
}
#endif
