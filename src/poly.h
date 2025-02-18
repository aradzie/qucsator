/*
 * poly.h - poly class definitions and implementation
 *
 * Copyright (C) 2005 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __POLY_H__
#define __POLY_H__

namespace qucs {

class poly
{
 public:
  poly (double _x, double _y)
    : x(_x), f0(_y) { f1 = f2 = 0; }
  poly (double _x, double _f0, double _f1)
    : x(_x), f0(_f0), f1(_f1) { f2 = 0; }
  poly (double _x, double _f0, double _f1, double _f2)
    : x(_x), f0(_f0), f1(_f1), f2(_f2) { }
  ~poly ()
    { }
  double eval (double _x) {
    double dx = _x - x; return f0 + dx * (f1 + dx * f2);
  }

  double x;  // x     - argument
  double f0; // f(x)  - value
  double f1; // f'(x) - 1st derivative
  double f2; // f"(x) - 2nd derivative
};

} // namespace qucs

#endif /* __POLY_H__ */
