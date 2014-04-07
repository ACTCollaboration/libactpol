#include <stdio.h>
#include <stdlib.h>

#include "actpol/actpol.h"
#include "actpol/dirfile.h"

int
main(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        const char *fn = argv[i];
        printf("%s\n", fn);

        ACTpolDirfile *dirfile = ACTpolDirfile_open(fn);

        int nsamples;
        double *ctime = ACTpolDirfile_read_double_channel(dirfile, "C_Time", &nsamples);
        printf("nsamples = %d\n", nsamples);
        printf("length = %g sec\n", ctime[nsamples-1]-ctime[0]);
    }
    return 0;
}
