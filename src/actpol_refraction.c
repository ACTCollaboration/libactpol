
#include <stdio.h>
#include <stdlib.h>

#include "actpol/actpol.h"
#include "actpol/util.h"

int
main(int argc, char *argv[])
{
    double freq_GHz = 230.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    printf("# temp[C] pressure[mbar] humidity[%] freq[GHz] alt[deg] ref[\"]\n");

    for (double temp_C = -20.; temp_C < 21.; temp_C += 10.)
    for (double pressure = 548.; pressure < 561.; pressure += 1.)
    for (double humidity = 0.; humidity < 1.1; humidity += 0.1)
    for (double alt_deg = 15.; alt_deg < 61.; alt_deg += 1.)
    {
        weather.temperature_C = temp_C;
        weather.pressure_mbar = pressure;
        weather.relative_humidity = humidity;

        double alt = deg2rad(alt_deg);
        double ref = actpol_refraction(&weather, freq_GHz, alt);
        double ref2 = actpol_refraction_alma366(&weather, freq_GHz, alt);
        double ref3 = actpol_refraction_ulich(&weather, freq_GHz, alt);
        double ref_arcsec = rad2arcsec(ref);

        printf("%g %g %g %g %g %g %g %g\n", temp_C,
            weather.pressure_mbar, weather.relative_humidity,
            freq_GHz, alt_deg, ref_arcsec, rad2arcsec(ref2), ref3);
    }

    return 0;
}

