
nobase_include_HEADERS = actpol/actpol.h \
                         actpol/array.h \
                         actpol/astro.h \
                         actpol/constants.h \
                         actpol/iers_bulletin_a.h \
                         actpol/map.h \
                         actpol/math.h \
                         actpol/oldact.h \
                         actpol/quaternion.h \
                         actpol/state.h \
                         actpol/vec3.h \
                         actpol/util.h

lib_LTLIBRARIES = libactpol.la
libactpol_la_SOURCES = array.c \
                       astro.c \
                       debug.h \
                       iers_bulletin_a.c \
                       map.c \
                       oldact.c \
                       quaternion.c \
                       state.c
libactpol_la_LIBADD = -lcfitsio -lslalib -lsofa_c -lwcs -lm

bin_PROGRAMS = actpol_refraction
actpol_refraction_LDADD = libactpol.la -lcfitsio -lslalib -lsofa_c -lwcs -lm

check_PROGRAMS = test
test_LDADD = libactpol.la -lcfitsio -lslalib -lsofa_c -lwcs -lm
TESTS = $(check_PROGRAMS)

if DIRFILE
nobase_include_HEADERS += actpol/dirfile.h actpol/getdata.h
libactpol_la_SOURCES += dirfile.c getdata.c
libactpol_la_LIBADD += -lzzip -lslim
actpol_refraction_LDADD += -lzzip -lslim
test_LDADD += -lzzip -lslim
endif

