/*
 * real.h - some real valued function definitions
 *
 * Copyright (C) 2008 Stefan Jahn <stefan@lkcc.org>
 * Copyright (C) 2014 Guilheme Brondani Torri <guitorri@gmail.com>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __REAL_H__
#define __REAL_H__

#include <cmath>

#include <constants.h>

/**It is preferred to add all used functions into the qucs namespace.
 * Doing so one is forced do think about compatibility instead of using std directly.
 * Inline is optional at this moment
 * \todo test if inline indeed performance improves (optimization flags should inline them anyway)
 */

namespace qucs {

//
// trigonometric
//
double    cos (const double);
double    sin (const double);
double    tan (const double);
double   acos (const double);
double   asin (const double);
double   atan (const double);
double  atan2 (const double, const double); //not used?, only for complex


//
// hyperbolic
//
double   cosh (const double);
double   sinh (const double);
double   tanh (const double);
double  acosh (const double); // c++11
double  asinh (const double); // c++11
double  atanh (const double); // c++11, not used?, only for complex


//
// exponential and logarithmic functions
//
double exp (const double);
double log (const double);
double log10 (const double);


//
// power functions
//
double pow (const double, const double );
double sqrt (const double );
double xhypot (const double, const double ); // same hypot in c++11?


//
// error functions
//
double erf(const double );


//
// rounding and remainder functions
//
double ceil(const double );
double floor(const double );
double fmod(const double ); //FIXME
double trunc(const double ); // c++11
double round(const double ); // c++11

//
// Qucs extra trigonometric helper
//
double coth (const double );
double sech (const double );
double cosech (const double );


//
// Qucs extra math functions
//
double  sqr (const double );
unsigned int sqr (unsigned int);
double  quadr (const double );

double rad2deg (const double );
double deg2rad (const double x );

/*!\brief Compute the third power of input */
static inline double cubic (const double x)  { return (x * x * x); }

/*!\brief Convert Celsius to Kelvin */
static inline double celsius2kelvin (const double x)  { return (x - K); }

/*!\brief Convert Kelvin to Celsius */
static inline double kelvin2celsius (const double x)  { return (x + K); }


//
// extra math functions
//
double limexp (const double);
double signum (const double);
double   sign (const double);
double   sinc (const double);
double    fix (const double);
double   step (const double);
unsigned int factorial (unsigned int);


//
// overload complex manipulations on reals
//
double   real (const double);
double   imag (const double);
double   norm (const double);
double   conj (const double);
double   abs (const double);

} // namespace qucs

#endif /* __REAL_H__ */
