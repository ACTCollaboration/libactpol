
#include <stdio.h>
#include <stdlib.h>

#include "actpol/actpol.h"

int
main(int argc, char *argv[])
{
    double freq_GHz = 150.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    printf("# temp[C] pressure[mbar] humidity[%] tropo-lapse[K/m] freq[GHz] alt[deg] ref[\"]\n");

    for (double temp_C = -20.; temp_C < 30; temp_C += 10.)
    for (double alt_deg = 35.; alt_deg < 61.; alt_deg += 1.)
    {
        weather.temperature_C = temp_C;
        double alt = deg2rad(alt_deg);
        double ref = actpol_refraction(&weather, freq_GHz, alt);
        double ref_arcsec = rad2arcsec(ref);

        printf("%g %g %g %g %g %g %g\n", temp_C,
            weather.pressure_mbar, weather.relative_humidity,
            weather.tropospheric_lapse_rate_K_per_m,
            freq_GHz, alt_deg, ref_arcsec);
    }

    return 0;
}
