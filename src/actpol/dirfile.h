
#pragma once

#include <stdbool.h>
#include <stdint.h>

struct FormatType;

typedef struct
{
    struct FormatType *format;
}
ACTpolDirfile;

ACTpolDirfile *
ACTpolDirfile_open(const char *filename);

void
ACTpolDirfile_close(ACTpolDirfile *dirfile);

bool
ACTpolDirfile_has_channel(const ACTpolDirfile *dirfile, const char *channel );

int16_t *
ACTpolDirfile_read_int16_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

uint16_t *
ACTpolDirfile_read_uint16_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

int32_t *
ACTpolDirfile_read_int32_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

uint32_t *
ACTpolDirfile_read_uint32_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

float *
ACTpolDirfile_read_float_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

double *
ACTpolDirfile_read_double_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples );

uint32_t
ACTpolDirfile_read_uint32_sample(const ACTpolDirfile *dirfile,
        const char *channelname, int index );

