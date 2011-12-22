
#include <assert.h>
#include <stdlib.h>

#include "actpol/array.h"
#include "actpol/astro.h"
#include "actpol/vec3.h"

ACTpolArray *
ACTpolArray_alloc(int ndets)
{
    ACTpolArray *array;
    array = (ACTpolArray *)malloc(sizeof(ACTpolArray));
    array->ndets = ndets;
    array->alt_offset = (double *)malloc(ndets*sizeof(double));
    array->az_offset = (double *)malloc(ndets*sizeof(double));
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
ACTpolArray_detector_alt_az(const ACTpolArray *array, int index,
        double boresight_alt, double boresight_az, double *alt, double *az)
{
    assert(index >= 0 && index < array->ndets);

    *alt = boresight_alt
         + array->boresight_offset_alt
         + array->alt_offset[index];
    *az = boresight_az
        + array->boresight_offset_az
        + array->az_offset[index];

    return 0;
}

int
ACTpolArray_compute_snapshot(ACTpolArray *array, Quaternion q, ACTpolArraySnapshot *snapshot)
{
    //double ref = refraction(&site, alt);
    double mat[3][3];
    Quaternion_to_matrix(q, mat);

    for (int i = 0; i != array->ndets; ++i) {
        double ref, rh[3], rc[3];
        actpol_ang2vec(-snapshot->az[i], snapshot->alt[i]-ref, rh);
        matrix_times_vec3(rc, mat, rh);
        actpol_vec2ang(rc, snapshot->ra+i, snapshot->dec+i);
    }

    return 0;
}

ACTpolArraySnapshot *
ACTpolArraySnapshot_alloc(ACTpolArray *array)
{
    ACTpolArraySnapshot *snapshot = (ACTpolArraySnapshot *)malloc(sizeof(ACTpolArraySnapshot));
    snapshot->array = array;
    snapshot->az = (double *)malloc(sizeof(double) * array->ndets);
    snapshot->alt = (double *)malloc(sizeof(double) * array->ndets);
    snapshot->ra = (double *)malloc(sizeof(double) * array->ndets);
    snapshot->dec = (double *)malloc(sizeof(double) * array->ndets);
    return snapshot;
}

void
ACTpolArraySnapshot_free(ACTpolArraySnapshot *snapshot)
{
    free(snapshot->dec);
    free(snapshot->ra);
    free(snapshot->alt);
    free(snapshot->az);
    free(snapshot);
}

