

nobase_include_HEADERS = actpol/actpol.h \
			 actpol/array.h \
			 actpol/astro.h \
			 actpol/constants.h \
			 actpol/dirfile.h \
		         actpol/getdata.h \
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
		      dirfile.c \
		      getdata.c \
		      iers_bulletin_a.c \
		      iers_bulletin_a.h \
		      map.c \
		      oldact.c \
		      quaternion.c \
		      state.c

bin_PROGRAMS = actpol_refraction
actpol_refraction_LDADD = libactpol.la -lcfitsio -lslalib -lsofa_c -lwcs -lzzip -lslim -lm

check_PROGRAMS = test
test_LDADD = libactpol.la -lcfitsio -lslalib -lsofa_c -lwcs -lzzip -lslim -lm
TESTS = $(check_PROGRAMS)

#EXTRA_DIST = iers_bulletin_a.dat \
#	     make_iers_bulletin_a_dot_c.py

