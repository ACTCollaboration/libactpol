//
// actpol/map.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

struct wcsprm;

typedef struct
{
    long naxis1, naxis2, npix;
    struct wcsprm *wcs;
    double *data;
}
ACTpolMap;

ACTpolMap *
ACTpolMap_alloc(long naxis1, long naxis2);

ACTpolMap *
ACTpolMap_new(double ra_deg_min, double ra_deg_max, double dec_deg_min, double dec_deg_max, double pixsize_deg);

ACTpolMap *
ACTpolMap_new_like(ACTpolMap *other_map);

ACTpolMap *
ACTpolMap_read_from_fits(const char *filename);

int
ACTpolMap_write_to_fits(ACTpolMap *map, const char *filename);

long
ACTpolMap_sky2pix(ACTpolMap *map, double ra_deg, double dec_deg);

long
ACTpolMap_sky2pix_cea_fast(ACTpolMap *map, double ra, double sindec);

int
ACTpolMap_pix2sky(ACTpolMap *map, long pix, double *ra_deg, double *dec_deg);

void
ACTpolMap_free(ACTpolMap *map);

#ifdef __cplusplus
}
#endif

