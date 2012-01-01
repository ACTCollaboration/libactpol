//
// actpol/oldact.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double latitude;
    double east_longitude;
    double elevation_m;
    double temperature_K;
    double pressure_mb;
    double relative_humidity;
}
ACTSite;

int
observed_altaz_to_mean_radec( const ACTSite *site, double freq_GHz,
        int n, const double ctime[], const double alt[], const double az[],
        double ra[], double dec[] );

#ifdef __cplusplus
}
#endif

