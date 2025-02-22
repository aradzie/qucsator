/*
 * sweep.cpp - variable sweep class implementation
 *
 * Copyright (C) 2004, 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include <cassert>
#include <cmath>
#include <cstring>

#include "object.h"
#include "sweep.h"
#include "vector.h"

namespace qucs {

sweep::sweep() : object() {
  type = SWEEP_UNKNOWN;
  data = nullptr;
  size = 0;
  text = nullptr;
  counter = 0;
}

sweep::sweep(const std::string &n) : object(n) {
  type = SWEEP_UNKNOWN;
  data = nullptr;
  size = 0;
  text = nullptr;
  counter = 0;
}

sweep::~sweep() {
  free(data);
  free(text);
}

// Returns the value at the given position.
double sweep::get(int idx) {
  assert(idx >= 0 && idx < size && data != nullptr);
  return data[idx];
}

// Sets the given value at the given position.
void sweep::set(int idx, double val) {
  assert(idx >= 0 && idx < size && data != nullptr);
  data[idx] = val;
}

/* Modifies the current size of the sweep.  If the new number of points
 * is larger than the current one it zeroes the new elements. */
void sweep::setSize(int points) {
  assert(points > 0);
  if (data != nullptr) {
    data = (double *)realloc(data, sizeof(double) * points);
    if (points > size) {
      memset(&data[size], 0, sizeof(double) * (points - size));
    }
  } else {
    data = (double *)malloc(sizeof(double) * points);
    memset(data, 0, sizeof(double) * points);
  }
  size = points;
  counter = 0;
}

char *sweep::toString() {
  free(text);
  if (data == nullptr || size == 0) {
    return (char *)"";
  }
  int len = 3 + size - 1;
  text = (char *)malloc(len);
  strcpy(text, "[");
  for (int i = 0; i < size; i++) {
    static char str[256]; // enough for a real number
    sprintf(str, "%g", get(i));
    text = (char *)realloc(text, len += strlen(str));
    strcat(text, str);
    if (i != size - 1) {
      strcat(text, ";");
    }
  }
  strcat(text, "]");
  return text;
}

/* Reverses the values order inside the sweep definition. */
void sweep::reverse() {
  if (data != nullptr && size > 0) {
    double *buf = (double *)malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++) {
      buf[i] = data[size - 1 - i];
    }
    free(data);
    data = buf;
  }
}

/* Returns the current sweep value and afterwards steps one value forward.
 * It wraps around at the end of sweep. */
double sweep::next() {
  double res = data[counter];
  if (++counter >= size) {
    counter = 0;
  }
  return res;
}

/* Returns the sweep value before the current value and thereby steps one value back.
 * It wraps around at the beginning of the sweep. */
double sweep::prev() {
  if (--counter < 0) {
    counter = size - 1;
  }
  return data[counter];
}

linsweep::linsweep() : sweep() { type = SWEEP_LINEAR; }

linsweep::linsweep(const std::string &n) : sweep(n) { type = SWEEP_LINEAR; }

/* Creates a linear stepped vector of values starting at the given start value,
 * ending with the given stop value and containing points elements. */
void linsweep::create(double start, double stop, int points) {
  vector v = linspace(start, stop, points);
  setSize(points);
  for (int i = 0; i < points; i++) {
    set(i, real(v.get(i)));
  }
}

linsweep::~linsweep() {}

logsweep::logsweep() : sweep() { type = SWEEP_LOGARITHMIC; }

logsweep::logsweep(const std::string &n) : sweep(n) { type = SWEEP_LOGARITHMIC; }

/* Creates a logarithmic stepped vector of values starting at the given start value,
 * ending with the given stop value and containing points elements. */
void logsweep::create(double start, double stop, int points) {
  vector v = logspace(start, stop, points);
  setSize(points);
  for (int i = 0; i < points; i++)
    set(i, real(v.get(i)));
}

logsweep::~logsweep() {}

consweep::consweep() : sweep() { type = SWEEP_CONSTANT; }

consweep::consweep(const std::string &n) : sweep(n) { type = SWEEP_CONSTANT; }

void consweep::create(double val) {
  setSize(1);
  set(0, val);
}

consweep::~consweep() {}

lstsweep::lstsweep() : sweep() { type = SWEEP_LIST; }

lstsweep::lstsweep(const std::string &n) : sweep(n) { type = SWEEP_LIST; }

/* Creates arbitrary values in a sweep containing points elements.
 * The actual values must be assigned using the set() function. */
void lstsweep::create(int points) { setSize(points); }

lstsweep::~lstsweep() {}

} // namespace qucs
