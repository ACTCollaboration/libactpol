
#include <assert.h>
#include <stdlib.h>

#include <slalib.h>

#include "actpol/math.h"
#include "actpol/oldact.h"
#include "iers_bulletin_a.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double secs2days( double s ) { return s/86400.; }
static inline double deg2rad( double deg ) { return deg*M_PI/180.; }
static inline double rad2deg( double rad ) { return rad*180./M_PI; }
static inline double arcsec2rad( double sec ) { return sec*M_PI/(180*3600); }

/*
 * Not strictly correct, because of leap seconds.
 */
double
convert_ctime_to_utc_mjd( double ctime )
{
    return secs2days(ctime) + 40587.0;
}

double
convert_utc_to_tt( double utc )
{
    return utc + secs2days(slaDtt(utc));
}

int
observed_altaz_to_mean_radec( const ACTSite *site, double freq_ghz,
        int n, const double ctime[],
        const double alt[], const double az[],
        double ra[], double dec[] )
{
    assert( n > 0 );
    assert( ctime != NULL );
    assert( alt != NULL );
    assert( az != NULL );
    assert( ra != NULL );
    assert( dec != NULL );

    int stat;
    double dut1, x, y;
    double amprms[21], aoprms[14];

    double utc = convert_ctime_to_utc_mjd( ctime[0] );

    stat = actpol_get_iers_bulletin_a( utc, &dut1, &x, &y );
    if ( stat != 0 )
        return stat;

    double wavelength_um = 299792.458/freq_ghz;

    slaAoppa( utc, dut1,
            site->east_longitude,
            site->latitude,
            site->elevation_m,
            arcsec2rad(x),
            arcsec2rad(y),
            site->temperature_K,
            site->pressure_mb,
            site->relative_humidity,
            wavelength_um,
            0.0065,     // tropospheric lapse rate [K/m]
            aoprms );

    double tt = convert_utc_to_tt( utc );
    // using tt instead of tdb
    slaMappa( 2000.0, tt, amprms );
    //printf( "freq = %g\n", freq_ghz );

    for ( int i = 0; i < n; i++ )
    {
        double observed_az = az[i];
        double observed_zenith = M_PI/2 - alt[i];
        double apparent_ra, apparent_dec;

        utc = convert_ctime_to_utc_mjd( ctime[i] );
        slaAoppat( utc, aoprms );

        slaOapqk( "A", observed_az, observed_zenith, aoprms,
                &apparent_ra, &apparent_dec );

        //printf( "apparent ra,dec = %g, %g\n", apparent_ra, apparent_dec );

        slaAmpqk( apparent_ra, apparent_dec, amprms, ra+i, dec+i );
    }

    return 0;
}

