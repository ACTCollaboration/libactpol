//
// actpol/math.h : libactpol header file
//
// 2012 Mike Nolta <mike@nolta.net>
//

#pragma once

#include <math.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
#define M_PI_2		1.57079632679489661923	// pi/2
#endif

#ifndef invsqrt
static double inline invsqrt( double x ) { return 1./sqrt(x); };
#endif

