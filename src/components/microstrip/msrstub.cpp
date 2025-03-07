/*
 * msrstub.cpp - microstrip radial stub class implementation
 *
 * Copyright (C) 2009 Stefan Jahn <stefan@lkcc.org>
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
#include "msrstub.h"

using namespace qucs;

msrstub::msrstub () : circuit (1) {
  type = CIR_MSRSTUB;
}

// Returns the microstrip radial stub reactance.
double msrstub::calcReactance (double r1, double r2,
				    double alpha, double er,
				    double h, double frequency) {

  double l0 = C0 / frequency;
  double W = (r1 + (r2 - r1) / 2) * deg2rad (alpha);
  double ereff = (er + 1.0) / 2 + (er - 1.0) /
    (2.0 * qucs::sqrt (1 + 10.0 * h / W));
  double k = 2.0 * pi * qucs::sqrt (ereff) / l0;
  double a = k * r1;
  double b = k * r2;
  double Z_0 = Z0 / qucs::sqrt (ereff) * qucs::sqrt (sqr (j0 (a)) + sqr (y0 (a))) /
    qucs::sqrt (sqr (j1 (a)) + sqr (y1 (a)));
  double theta_1 = qucs::atan (y0 (a) / j0 (a));
  //  double theta_2 = atan (y0 (b) / j0 (b));
  double phi_1 = qucs::atan (-j1 (a) / y1 (a));
  double phi_2 = qucs::atan (-j1 (b) / y1 (b));

  double X1 = h * Z_0 / (2.0 * pi * r1) * 360.0 / alpha *
    qucs::cos (theta_1 - phi_2) / qucs::sin (phi_1 - phi_2);

  return X1;
}

void msrstub::calcSP (double frequency) {
  setS (NODE_1, NODE_1, ztor (calcZ (frequency)));
}

nr_complex_t msrstub::calcZ (double frequency) {

  /* get properties of this component */
  double r1 = getPropertyDouble ("ri");
  double r2 = getPropertyDouble ("ro");
  double al = getPropertyDouble ("alpha");

  /* get properties of the substrate */
  substrate * subst = getSubstrate ();
  double er    = subst->getPropertyDouble ("er");
  double h     = subst->getPropertyDouble ("h");

  return nr_complex_t (0, calcReactance (r1, r2, al, er, h, frequency));
}

void msrstub::initDC (void) {
  allocMatrixMNA ();
  setY (NODE_1, NODE_1, 0);
}

void msrstub::calcAC (double frequency) {
  setY (NODE_1, NODE_1, 1.0 / calcZ (frequency));
}

// properties
PROP_REQ [] = {
  { "ri", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "ro", PROP_REAL, { 10e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "alpha", PROP_REAL, { 90, PROP_NO_STR }, PROP_RNGII (0, 180) },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t msrstub::cirdef =
  { "MRSTUB", 1, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
