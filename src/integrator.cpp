/*
 * integrator.cpp - integrator class implementation
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#include "integrator.h"

namespace qucs {

integrator::integrator()
    : mode(0), order(0), coefficients(nullptr), integrate_func(nullptr), conductor_func(nullptr) {}

/* Evaluates the state of the integration-using component
 * and runs the appropriate integrator function. */
void integrator::integrate(int qstate, double cap, double &geq, double &ceq) {
  int cstate = qstate + 1;
  if (mode & MODE_INIT) {
    fillState(qstate, getState(qstate));
  }
  (*integrate_func)(this, qstate, cap, geq, ceq);
  if (mode & MODE_INIT) {
    fillState(cstate, getState(cstate));
  }
}

/* Runs the appropriate conductor function. */
void integrator::conductor(double cap, double &geq) { (*conductor_func)(this, cap, geq); }

} // namespace qucs
