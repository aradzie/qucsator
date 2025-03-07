/*
 * vac.cpp - AC voltage source class implementation
 *
 * Copyright (C) 2003, 2004, 2006, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "vac.h"

using namespace qucs;

vac::vac () : circuit (2) {
  type = CIR_VAC;
  setVSource (true);
  setVoltageSources (1);
}

void vac::initSP (void) {
  allocMatrixS ();
  setS (NODE_1, NODE_1, 0.0);
  setS (NODE_1, NODE_2, 1.0);
  setS (NODE_2, NODE_1, 1.0);
  setS (NODE_2, NODE_2, 0.0);
}

void vac::initDC (void) {
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void vac::initAC (void) {
  initDC ();
  double a = getPropertyDouble ("U");
  double p = getPropertyDouble ("Phase");
  setE (VSRC_1, qucs::polar (a, deg2rad (p)));
}

void vac::initTR (void) {
  initDC ();
}

void vac::calcTR (double t) {
  double f = getPropertyDouble ("f");
  double p = getPropertyDouble ("Phase");
  double d = getPropertyDouble ("Theta");
  double a = getPropertyDouble ("U");
  double s = getNet()->getSrcFactor ();
  double o = 2 * pi * f;
  double T = p / f / 360;
  double u = s * a * std::exp (-(t + T) * d * f) * std::sin (o * t + deg2rad (p));
  setE (VSRC_1, u);
}

void vac::initHB (void) {
  setVoltageSources (1);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void vac::calcHB (double frequency) {
  double f = getPropertyDouble ("f");
  if (f == frequency) {
    double a = getPropertyDouble ("U");
    double p = getPropertyDouble ("Phase");
    setE (VSRC_1, qucs::polar (a, deg2rad (p)));
  }
  else {
    setE (VSRC_1, 0);
  }
}

// properties
PROP_REQ [] = {
  { "U", PROP_REAL, { 1, PROP_NO_STR }, PROP_NO_RANGE }, PROP_NO_PROP };
PROP_OPT [] = {
  { "Phase", PROP_REAL, { 0, PROP_NO_STR }, PROP_RNGII (-360, 360) },
  { "Theta", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "f", PROP_REAL, { 1e9, PROP_NO_STR }, PROP_POS_RANGE },
  PROP_NO_PROP };
struct define_t vac::cirdef =
  { "Vac", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
