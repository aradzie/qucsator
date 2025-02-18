/*
 * msvia.cpp - microstrip via hole class implementation
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

#include "component.h"
#include "substrate.h"
#include "msvia.h"

using namespace qucs;

msvia::msvia () : circuit (2) {
  R = 0;
  Z = 0;
  type = CIR_MSVIA;
}

void msvia::calcNoiseSP (double) {
  // calculate noise correlation matrix
  double T = getPropertyDouble ("Temp");
  double f = celsius2kelvin (T) * 4.0 * real (Z) * z0 / norm (4.0 * z0 + Z) / T0;
  setN (NODE_1, NODE_1, +f); setN (NODE_2, NODE_2, +f);
  setN (NODE_1, NODE_2, -f); setN (NODE_2, NODE_1, -f);
}

void msvia::initSP (void) {
  allocMatrixS ();
  R = calcResistance ();
}

void msvia::calcSP (double frequency) {
  // calculate s-parameters
  Z = calcImpedance (frequency);
  nr_complex_t z = Z / z0;
  setS (NODE_1, NODE_1, z / (z + 2.0));
  setS (NODE_2, NODE_2, z / (z + 2.0));
  setS (NODE_1, NODE_2, 2.0 / (z + 2.0));
  setS (NODE_2, NODE_1, 2.0 / (z + 2.0));
}

nr_complex_t msvia::calcImpedance (double frequency) {
  // fetch substrate and component properties
  substrate * subst = getSubstrate ();
  double h   = subst->getPropertyDouble ("h");
  double t   = subst->getPropertyDouble ("t");
  double rho = subst->getPropertyDouble ("rho");
  double r   = getPropertyDouble ("D") / 2;

  // check frequency validity
  if (frequency * h >= 0.03 * C0) {
    logprint (LOG_ERROR, "WARNING: Model for microstrip via hole defined for "
	      "freq/C0*h < 0.03 (is %g)\n", frequency / C0 * h);
  }

  // create Z-parameter
  double fs  = pi * MU0 * sqr (t) / rho;
  double res = R * std::sqrt (1 + frequency * fs);
  double a   = std::sqrt (sqr (r) + sqr (h));
  double ind = MU0 * (h * std::log ((h + a) / r) + 1.5 * (r - a));
  return Z = nr_complex_t (res, frequency * ind);
}

double msvia::calcResistance (void) {
  // fetch substrate and component properties
  substrate * subst = getSubstrate ();
  double h   = subst->getPropertyDouble ("h");
  double t   = subst->getPropertyDouble ("t");
  double rho = subst->getPropertyDouble ("rho");
  double r   = getPropertyDouble ("D") / 2;
  double v   = h / pi / (sqr (r) - sqr (r - t));
  return R = rho * v;
}

void msvia::initDC (void) {
  double r = calcResistance ();

  // for non-zero resistances usual MNA entries
  if (r != 0.0) {
    double g = 1.0 / r;
    setVoltageSources (0);
    allocMatrixMNA ();
    setY (NODE_1, NODE_1, +g); setY (NODE_2, NODE_2, +g);
    setY (NODE_1, NODE_2, -g); setY (NODE_2, NODE_1, -g);
  }
  // for zero resistances create a zero voltage source
  else {
    setVoltageSources (1);
    setInternalVoltageSource (1);
    allocMatrixMNA ();
    clearY ();
    voltageSource (VSRC_1, NODE_1, NODE_2);
  }
}

void msvia::initAC (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
  R = calcResistance ();
}

void msvia::calcAC (double frequency) {
  nr_complex_t y = 1.0 / calcImpedance (frequency);
  setY (NODE_1, NODE_1, +y); setY (NODE_2, NODE_2, +y);
  setY (NODE_1, NODE_2, -y); setY (NODE_2, NODE_1, -y);
}

void msvia::calcNoiseAC (double) {
  // calculate noise current correlation matrix
  double y = real (1.0 / Z);
  double T = getPropertyDouble ("Temp");
  double f = celsius2kelvin (T) / T0 * 4.0 * y;
  setN (NODE_1, NODE_1, +f); setN (NODE_2, NODE_2, +f);
  setN (NODE_1, NODE_2, -f); setN (NODE_2, NODE_1, -f);
}

// properties
PROP_REQ [] = {
  { "D", PROP_REAL, { 100e-6, PROP_NO_STR }, PROP_POS_RANGE },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  PROP_NO_PROP };
struct define_t msvia::cirdef =
  { "MVIA", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
