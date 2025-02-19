/*
 * spline.h - spline class definitions
 *
 * Copyright (C) 2005, 2006, 2009 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __SPLINE_H__
#define __SPLINE_H__

#include <vector>

#include "tvector.h"
#include "vector.h"

namespace qucs {

// Types of boundary conditions.
enum spline_boundary_type {
  SPLINE_BC_UNKNOWN = -1,
  SPLINE_BC_NATURAL, // natural splines -- zero derivatives
  SPLINE_BC_CLAMPED, // endpoint derivatives given
  SPLINE_BC_PERIODIC // periodic splines
};

class vector;
class poly;

class spline {
public:
  spline();
  spline(int);
  spline(tvector<double>, tvector<double>);
  spline(qucs::vector, qucs::vector);
  spline(std::vector<double>, std::vector<double>);
  ~spline();

  void vectors(qucs::vector, qucs::vector);
  void vectors(tvector<double>, tvector<double>);
  void vectors(std::vector<double>, std::vector<double>);
  void vectors(double *, double *, int);
  void construct(void);
  poly evaluate(double);
  void setBoundary(int b) { boundary = b; }
  void setDerivatives(double l, double r) {
    d0 = l;
    dn = r;
  }

private:
  double *upper_bound(double *, double *, double);
  void realloc(int);

private:
  double *x;
  double *f0;
  double *f1;
  double *f2;
  double *f3;
  double d0, dn;
  int n;
  int boundary;
};

} // namespace qucs

#endif /* __SPLINE_H__ */
