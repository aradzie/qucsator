/*
 * states.cpp - save-state variable class implementation
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
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

#include <cstdlib>

#include "states.h"

namespace qucs {

constexpr int STATE_SHIFT = 3;
constexpr int STATE_NUM = 8;
constexpr int STATE_MASK = 7;

states::states() {
  size = 0;
  values = nullptr;
  pos = 0;
}

states::~states() {
  if (values) {
    free(values);
  }
}

void states::setStates(int const n) {
  size = n;
  if (values) {
    free(values);
  }
  values = (double *)calloc(size, sizeof(double) * STATE_NUM);
  pos = 0;
}

int states::getStates() const { return size; }

/* Applies the given value to a save-state variable through all history values. */
void states::fillState(const int stateId, double value) {
  for (int i = 0; i < STATE_NUM; i++) {
    values[(stateId << STATE_SHIFT) + i] = value;
  }
}

/* Applies the given value to a save-state variable.
   Higher positions mean earlier states.
   By default the function sets the current state of the save-state variable. */
void states::setState(const int stateId, const double value, const int offset) {
  const int index = (offset + pos) & STATE_MASK;
  values[(stateId << STATE_SHIFT) + index] = value;
}

/* Returns a save-state variable at the given position.
   Higher positions mean earlier states.
   By default the function returns the current state of the save-state variable. */
double states::getState(const int stateId, const int offset) const {
  const int index = (offset + pos) & STATE_MASK;
  return values[(stateId << STATE_SHIFT) + index];
}

// Shifts one state forward.
void states::nextState() {
  if (pos > 0) {
    pos--;
  } else {
    pos = STATE_NUM - 1;
  }
}

} // namespace qucs
