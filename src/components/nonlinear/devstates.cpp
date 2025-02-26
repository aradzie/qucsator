/*
 * devstates.cpp - device state class implementation
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

#include <cstdlib>

#include "devstates.h"

namespace qucs {

devstates::devstates() {
  nstates = 0;
  nvars = 0;
  states = nullptr;
  nstate = 0;
  pstate = nullptr;
}

devstates::devstates(int vars, int states) { deviceStates(vars, states); }

devstates::~devstates() { free(states); }

/* Initializes the device state class instance containing the
   specified number of variables and states. */
void devstates::deviceStates(int vars, int stats) {
  nvars = vars;
  nstates = stats;
  free(states);
  states = (double *)malloc(sizeof(double) * nvars * nstates);
  nstate = 0;
  pstate = states;
}

// Returns the number of states.
int devstates::deviceStates() { return nstates; }

// Sets the current state.
void devstates::deviceState(int state) {
  nstate = state;
  pstate = &states[nvars * nstate];
}

// Returns the current state.
int devstates::deviceState() { return nstate; }

// Access operator for the given variable in the current state.
double devstates::operator()(int var) const { return pstate[var]; }

// Reference access operator for the given variable in the current state.
double &devstates::operator()(int var) { return pstate[var]; }

// Returns the given variable in the current state.
double devstates::deviceVar(int var) const { return pstate[var]; }

// Returns a reference to the given variable in the current state.
double &devstates::deviceVar(int var) { return pstate[var]; }

} // namespace qucs
