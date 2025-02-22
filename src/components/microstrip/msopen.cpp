/*
 * msopen.cpp - microstrip open end class implementation
 *
 * Copyright (C) 2004, 2008 Stefan Jahn <stefan@lkcc.org>
 * Copyright (C) 2004 Michael Margraf <Michael.Margraf@alumni.TU-Berlin.DE>
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
#include "msline.h"
#include "msopen.h"

using namespace qucs;

msopen::msopen () : circuit (1) {
  type = CIR_MSOPEN;
}

// Returns the microstrip open end capacitance.
double msopen::calcCend (double frequency, double W,
			      double h, double t, double er,
			      const char * const SModel, const char * const DModel,
			      const char * const Model) {

  double ZlEff, ErEff, WEff, ZlEffFreq, ErEffFreq;
  msline::analyseQuasiStatic (W, h, t, er, SModel, ZlEff, ErEff, WEff);
  msline::analyseDispersion  (WEff, h, er, ZlEff, ErEff, frequency, DModel,
			      ZlEffFreq, ErEffFreq);

  W /= h;
  double dl = 0;
  /* Kirschning, Jansen and Koster */
  if (!strcmp (Model, "Kirschning")) {
    double Q6 = qucs::pow (ErEffFreq, 0.81);
    double Q7 = qucs::pow (W, 0.8544);
    double Q1 = 0.434907 *
      (Q6 + 0.26) / (Q6 - 0.189) * (Q7 + 0.236) / (Q7 + 0.87);
    double Q2 = qucs::pow (W, 0.371) / (2.358 * er + 1.0) + 1.0;
    double Q3 = qucs::atan (0.084 * qucs::pow (W, 1.9413 / Q2)) *
      0.5274 / qucs::pow (ErEffFreq, 0.9236) + 1.0;
    double Q4 = 0.0377 * (6.0 - 5.0 * qucs::exp (0.036 * (1.0 - er))) *
      qucs::atan (0.067 * qucs::pow (W, 1.456)) + 1.0;
    double Q5 = 1.0 - 0.218 * qucs::exp (-7.5 * W);
    dl = Q1 * Q3 * Q5 / Q4;
  }
  /* Hammerstad */
  else if (!strcmp (Model, "Hammerstad")) {
    dl = 0.102 * (W + 0.106) / (W + 0.264) *
      (1.166 + (er + 1) / er * (0.9 + qucs::log (W + 2.475)));
  }
  return dl * h * qucs::sqrt (ErEffFreq) / C0 / ZlEffFreq;
}

void msopen::calcSP (double frequency) {
  setS (NODE_1, NODE_1, ztor (1.0 / calcY (frequency)));
}

nr_complex_t msopen::calcY (double frequency) {

  /* how to get properties of this component, e.g. W */
  double W = getPropertyDouble ("W");
  const char * SModel = getPropertyString ("MSModel");
  const char * DModel = getPropertyString ("MSDispModel");
  const char * Model  = getPropertyString ("Model");

  /* how to get properties of the substrate, e.g. Er, H */
  substrate * subst = getSubstrate ();
  double er    = subst->getPropertyDouble ("er");
  double h     = subst->getPropertyDouble ("h");
  double t     = subst->getPropertyDouble ("t");

  /* local variables */
  nr_complex_t y;
  double o = 2 * pi * frequency;

  /* Alexopoulos and Wu */
  if (!strcmp (Model, "Alexopoulos")) {
    double ZlEff, ErEff, WEff, ZlEffFreq, ErEffFreq;
    msline::analyseQuasiStatic (W, h, t, er, SModel, ZlEff, ErEff, WEff);
    msline::analyseDispersion  (WEff, h, er, ZlEff, ErEff, frequency, DModel,
				ZlEffFreq, ErEffFreq);

    if (fabs (er - 9.9) > 0.2) {
      logprint (LOG_ERROR, "WARNING: Model for microstrip open end defined "
		"for er = 9.9 (er = %g)\n", er);
    }

    double c1, c2, l2, r2;
    c1 = (1.125 * qucs::tanh (1.358 * W / h) - 0.315) *
      h / 2.54e-5 / 25 / ZlEffFreq * 1e-12;
    c2 = (6.832 * qucs::tanh (0.0109 * W / h) + 0.919) *
      h / 2.54e-5 / 25 / ZlEffFreq * 1e-12;
    l2 = (0.008285 * qucs::tanh (0.5665 * W / h) + 0.0103) *
      h / 2.54e-5 / 25 * ZlEffFreq * 1e-9;
    r2 = (1.024 * qucs::tanh (2.025 * W / h)) * ZlEffFreq;
    y = nr_complex_t (0, c1 * o) + 1.0 / nr_complex_t (r2, l2 * o - 1 / c2 / o);
  }
  else {
    double c = calcCend (frequency, W, h, t, er, SModel, DModel, Model);
    y = nr_complex_t (0, c * o);
  }
  return y;
}

void msopen::initDC (void) {
  allocMatrixMNA ();
  setY (NODE_1, NODE_1, 0);
}

void msopen::calcAC (double frequency) {
  setY (NODE_1, NODE_1, calcY (frequency));
}

// properties
PROP_REQ [] = {
  { "W", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  { "MSDispModel", PROP_STR, { PROP_NO_VAL, "Kirschning" }, PROP_RNG_DIS },
  { "MSModel", PROP_STR, { PROP_NO_VAL, "Hammerstad" }, PROP_RNG_MOD },
  { "Model", PROP_STR, { PROP_NO_VAL, "Kirschning" },
    PROP_RNG_STR3 ("Kirschning", "Hammerstad", "Alexopoulos") },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t msopen::cirdef =
  { "MOPEN", 1, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
