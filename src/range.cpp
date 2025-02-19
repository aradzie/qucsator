/*
 * range.cpp - range class implementation
 *
 * Copyright (C) 2006 Stefan Jahn <stefan@lkcc.org>
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

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "range.h"

namespace qucs {

range::range() {
  il = ih = '.';
  l = h = 0.0;
  txt = nullptr;
}

range::range(char ilo, double lo, double hi, char ihi) {
  il = ilo;
  ih = ihi;
  if (lo > hi) {
    h = lo;
    l = hi;
  } else {
    l = lo;
    h = hi;
  }
  txt = nullptr;
}

range::range(const range &r) {
  txt = r.txt ? strdup(r.txt) : nullptr;
  il = r.il;
  ih = r.ih;
  l = r.l;
  h = r.h;
}

/* Checks whether the given value is outside the range. */
bool range::outside(double value) { return !inside(value); }

/* Checks whether the given value is inside the range. */
bool range::inside(double value) {
  int err = 0;
  if (il == '[' && (value < l))
    err++;
  if (il == ']' && !(value > l))
    err++;
  if (ih == '[' && !(value < h))
    err++;
  if (ih == ']' && (value > h))
    err++;
  return err == 0;
}

range::~range() { free(txt); }

char *range::toString() {
  char str[64];
  sprintf(str, "%c%g,%g%c", il, l, h, ih);
  free(txt);
  txt = strdup(str);
  return txt;
}

} // namespace qucs
