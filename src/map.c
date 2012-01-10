
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fitsio.h>
#include <wcslib/wcslib.h>

#include "actpol/map.h"

static struct wcsprm *
alloc_wcs()
{
    wcserr_enable(1);
    struct wcsprm *wcs = (struct wcsprm *) malloc(sizeof(struct wcsprm));
    wcs->flag = -1;
    return wcs;
}

static struct wcsprm *
new_wcs(int naxis)
{
    struct wcsprm *wcs = alloc_wcs();
    int stat = wcsini(1, naxis, wcs);
    assert(stat == 0);
    return wcs;
}

static struct wcsprm *
copy_wcs(const struct wcsprm *iwcs)
{
    struct wcsprm *wcs = alloc_wcs();
    int stat = wcssub(1, iwcs, 0, 0, wcs);
    assert(stat == 0);
    return wcs;
}

static void
free_wcs(struct wcsprm *wcs)
{
    int stat = wcsfree(wcs);
    assert(stat == 0);
    free(wcs);
}

ACTpolMap *
ACTpolMap_alloc(long naxis1, long naxis2)
{
    ACTpolMap *map = (ACTpolMap *) malloc(sizeof(ACTpolMap));
    map->naxis1 = naxis1;
    map->naxis2 = naxis2;
    map->npix = map->naxis1*map->naxis2;
    map->wcs = NULL;
    map->data = (double *) malloc(sizeof(double)*map->npix);
    return map;
}

ACTpolMap *
ACTpolMap_new(double ra_min, double ra_max, double dec_min, double dec_max, double pixsize)
{
    long naxis1 = (long) round((ra_max - ra_min)/pixsize);
    long naxis2 = (long) round((dec_max - dec_min)/pixsize);
    assert(naxis1 > 0);
    assert(naxis2 > 0);

    struct wcsprm *wcs = new_wcs(2);

    strcpy(wcs->ctype[0], "RA---CEA");
    strcpy(wcs->ctype[1], "DEC--CEA");
    wcs->crval[0] = 0.5*(ra_min + ra_max);
    wcs->crval[1] = 0.5*(dec_min + dec_max);
    wcs->cdelt[0] = -pixsize;
    wcs->cdelt[1] = pixsize;
    wcs->crpix[0] = 0.5*naxis1;
    wcs->crpix[1] = 0.5*naxis2;
    wcs->pc[0] = 1.;
    wcs->pc[1] = 0.;
    wcs->pc[2] = 0.;
    wcs->pc[3] = 1.;
    wcs->npv = 1;
    wcs->pv[0].i = 2;
    wcs->pv[0].m = 1;
    wcs->pv[0].value = 1.;

    if (wcsset(wcs)) {
        wcsperr(wcs, "libactpol");
        free_wcs(wcs);
        return NULL;
    }

    ACTpolMap *map = ACTpolMap_alloc(naxis1, naxis2);
    map->wcs = wcs;
    return map;
}

ACTpolMap *
ACTpolMap_new_like(ACTpolMap *other_map)
{
    ACTpolMap *map = ACTpolMap_alloc(other_map->naxis1, other_map->naxis2);
    map->wcs = copy_wcs(other_map->wcs);
    return map;
}

long
ACTpolMap_sky2pix(ACTpolMap *map, double ra, double dec)
{
    double world[2] = {ra, dec};
    double phi[2], theta[2], imgcrd[2], pixcrd[2];
    int wstat[2] = {0, 0};
    wcss2p(map->wcs, 1, 2, world, phi, theta, imgcrd, pixcrd, wstat);
    if (wstat[0] || wstat[1]) {
        wcsperr(map->wcs, "libactpol");
        return -1;
    }

    //printf("sky2pix pixcrd %g %g\n", pixcrd[0], pixcrd[1]);
    long pix1 = (long) round(pixcrd[0] - 1.);
    long pix2 = (long) round(pixcrd[1] - 1.);
    if (pix1 >= 0 && pix1 < map->naxis1 && pix2 >= 0 && pix2 < map->naxis2)
        return pix1 + map->naxis1*pix2;
    return -1;
}

long
ACTpolMap_sky2pix_cea_fast(ACTpolMap *map, double ra, double sindec)
{
    double world[2] = {ra, sindec};
    double imgcrd[2], pixcrd[2];
    long pix[2];
    for (int i = 0; i < 2; i++)
    {
        imgcrd[i] = fmod(world[i]*180./M_PI - map->wcs->crval[i], 360.);
        //imgcrd[i] = fmod(world[i] - map->wcs->crval[i], 360.);
        if (imgcrd[i] > 180.)
            imgcrd[i] -= 360.;
        else if (imgcrd[i] < -180.)
            imgcrd[i] += 360.;
        pixcrd[i] = imgcrd[i]/map->wcs->cdelt[i] + map->wcs->crpix[i];
        pix[i] = (long) round(pixcrd[i] - 1.);
    }

    //printf("fast: %g,%g %g,%g %g,%g\n", world[0], world[1], imgcrd[0], imgcrd[1], pixcrd[0], pixcrd[1]);

    if (pix[0] >= 0 && pix[0] < map->naxis1 && pix[1] >= 0 && pix[1] < map->naxis2)
        return pix[0] + map->naxis1*pix[1];
    return -1;
}

int
ACTpolMap_pix2sky(ACTpolMap *map, long pix, double *ra, double *dec)
{
    assert(pix >= 0 && pix < map->npix);

    double imgcrd[2], phi[2], theta[2], world[2];
    long pix1 = pix % map->naxis1;
    long pix2 = (pix - pix1)/map->naxis1;
    double pixcrd[2] = {1. + pix1, 1. + pix2};
    //printf("pix2sky pixcrd %g %g\n", pixcrd[0], pixcrd[1]);
    int wstat[2] = {0, 0};
    wcsp2s(map->wcs, 1, 2, pixcrd, imgcrd, phi, theta, world, wstat);
    if (wstat[0] || wstat[1]) {
        //printf("err! wstat = %d %d\n", wstat[0], wstat[1]);
        wcsperr(map->wcs, "libactpol");
        return -1;
    }

    *ra = world[0];
    *dec = world[1];
    return 0;
}

int
ACTpolMap_write_to_fits(ACTpolMap *map, const char *filename)
{
    assert(map != NULL);
    assert(filename != NULL);
    assert(map->npix == map->naxis1*map->naxis2);

    int fstat = 0; // must be initialized
    fitsfile *fp = NULL;
    fits_create_file(&fp, filename, &fstat);
    if (fstat) {
        fits_report_error(stderr, fstat);
        return fstat;
    }

    long naxes[2] = {map->naxis1, map->naxis2};
    fits_create_img(fp, -64, 2, naxes, &fstat);

    char *header;
    int nkeyrec;
    wcshdo(0, map->wcs, &nkeyrec, &header);
    for (int i = 0; i != nkeyrec; ++i)
        fits_write_record(fp, header+i*80, &fstat);
    free(header);

    long fpixel[2] = {1, 1};
    fits_write_pix(fp, TDOUBLE, fpixel, map->npix, map->data, &fstat);

    fits_close_file(fp, &fstat);

    if (fstat)
        fits_report_error(stderr, fstat);
    return fstat;
}

ACTpolMap *
ACTpolMap_read_from_fits(const char *filename)
{
    int fstat = 0; // must be initialized
    fitsfile *fp = NULL;
    fits_open_file(&fp, filename, READONLY, &fstat);
    if (fstat) {
        fits_report_error(stderr, fstat);
        return NULL;
    }

    int naxis;
    fits_get_img_dim(fp, &naxis, &fstat);
    assert(naxis == 2);
    long naxes[2];
    fits_get_img_size(fp, 2, naxes, &fstat);
    if (fstat) {
        fits_report_error(stderr, fstat);
        return NULL;
    }

    char *header;
    int nkeys, nreject, nwcs;
    struct wcsprm *wcs;
    fits_hdr2str(fp, 1, NULL, 0, &header, &nkeys, &fstat);
    int wstat = wcspih(header, nkeys, 0, 2, &nreject, &nwcs, &wcs);
    free(header);
    if (wstat || nreject || nwcs != 1) {
        return NULL;
    }

    ACTpolMap *map = ACTpolMap_alloc(naxes[0], naxes[1]);
    map->wcs = wcs;

    long fpixel[2] = {1, 1};
    int anynull;
    fits_read_pix(fp, TDOUBLE, fpixel, map->npix, NULL, map->data, &anynull, &fstat);
    fits_close_file(fp, &fstat);
    if (fstat)
        fits_report_error(stderr, fstat);

    return map;
}

void
ACTpolMap_free(ACTpolMap *map)
{
    free_wcs(map->wcs);
    free(map->data);
    free(map);
}

