
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
			 actpol/vec3.h

lib_LIBRARIES = libactpol.a
libactpol_a_SOURCES = array.c \
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
actpol_refraction_LDADD = libactpol.a -lcfitsio -lslalib -lsofa_c -lwcs -lzzip -lslim

check_PROGRAMS = test
test_LDADD = libactpol.a -lcfitsio -lslalib -lsofa_c -lwcs -lzzip -lslim
TESTS = $(check_PROGRAMS)

#EXTRA_DIST = iers_bulletin_a.dat \
#	     make_iers_bulletin_a_dot_c.py

