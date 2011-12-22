
nobase_include_HEADERS = actpol/actpol.h \
			 actpol/array.h \
			 actpol/astro.h \
			 actpol/constants.h \
			 actpol/map.h \
			 actpol/quaternion.h

lib_LIBRARIES = libactpol.a
libactpol_a_SOURCES = array.c \
		      astro.c \
		      iers_bulletin_a.c \
		      iers_bulletin_a.h \
		      map.c \
		      quaternion.c

check_PROGRAMS = test
test_LDADD = libactpol.a -lcfitsio -lslalib -lsofa_c -lwcs
#test_map_LDADD = libactpol.a -lcfitsio -lwcs
#test_quaternion_LDADD = libactpol.a -lsofa_c
#test_astro_LDADD = libactpol.a -lsofa_c -lslalib
TESTS = $(check_PROGRAMS)

#EXTRA_DIST = iers_bulletin_a.dat \
#	     make_iers_bulletin_a_dot_c.py

