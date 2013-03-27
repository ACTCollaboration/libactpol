
libactpol
=========

Usage
-----

### Setup an array

    #include <actpol/actpol.h>

    ACTpolArray *array = ACTpolArray_alloc(nhorns);
    ACTpolArray_init(array, freq_GHz, array_center_x, array_center_y);

    for (int i = 0; i < nhorns; i++) {
        ...
        ACTpolFeedhorn_init(array->horn+i, x, y, a);
    }

### Compute array coordinates

    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_RA_SINDEC);

    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);

    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    for (int i = 0; i < nscans; i++)
    {
      ACTpolScan scan;
      ACTpolScan_init(&scan, scan_alt, scan_az, scan_throw);

      // only needs to be called when scan or weather changes
      ACTpolArrayCoords_update_refraction(coords, &scan, &weather);
      ...
      for (j = 0; j < nsamples; j++)
      {
        ...
        ACTpolState_update(state, unixtime, boresight_alt, boresight_az);
        ACTpolArrayCoords_update(coords, state);

        for (k = 0; k < nhorns; k++)
        {
          ACTpolFeedhornCoords *fc = coords->horn + k;
          printf("%g %g\n", fc->a, fc->b); // a = ra, b = sin(dec)
          printf("%g %g\n", fc->sin2gamma1, fc->cos2gamma1);
          printf("%g %g\n", fc->sin2gamma2, fc->cos2gamma2);
        }
      }
    }

### Read dirfile

    #include <actpol/dirfile.h>

    ACTpolDirfile *dirfile = ACTpolDirfile_open("/path/to/dirfile");

    if (ACTpolDirfile_has_channel(dirfile, "Enc_Az_Deg")) {
      int nsamples;
      double *az = ACTpolDirfile_read_double_channel(dirfile, "Enc_Az_Deg", &nsamples);
      ...
    }

    ACTpolDirfile_close(dirfile);

### Create and write map

    ACTpolMap *map = ACTpolMap_new(ra_min, ra_max, dec_min, dec_max, pixsize);
    for (...)
    {
      ...
      //long pix = ACTpolMap_sky2pix(map, ra, dec);
      long pix = ACTpolMap_sky2pix_cea_fast(map, ra, sindec);
      map->data[pix] += ...
    }
    ACTpolMap_write_to_fits(map, "filename");
    ACTpolMap_free(map);

