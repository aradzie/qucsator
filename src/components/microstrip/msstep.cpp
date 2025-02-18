/*
 * msstep.cpp - microstrip impedance step class implementation
 *
 * Copyright (C) 2004, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "msstep.h"

using namespace qucs;

msstep::msstep () : circuit (2) {
  type = CIR_MSSTEP;
}

void msstep::calcSP (double frequency) {
  setMatrixS (ztos (calcMatrixZ (frequency)));
}

matrix msstep::calcMatrixZ (double frequency) {

  /* how to get properties of this component, e.g. W */
  double W1 = getPropertyDouble ("W1");
  double W2 = getPropertyDouble ("W2");
  const char * SModel = getPropertyString ("MSModel");
  const char * DModel = getPropertyString ("MSDispModel");

  /* how to get properties of the substrate, e.g. Er, H */
  substrate * subst = getSubstrate ();
  double er    = subst->getPropertyDouble ("er");
  double h     = subst->getPropertyDouble ("h");
  double t     = subst->getPropertyDouble ("t");

  // compute parallel capacitance
  double t1 = std::log10 (er);
  double t2 = W1 / W2;
  double Cs = std::sqrt (W1 * W2) *
    (t2 * (10.1 * t1 + 2.33) - 12.6 * t1 - 3.17);

  // compute series inductance
  t1 = std::log10 (t2);
  t2 = t2 - 1;
  double Ls = h * (t2 * (40.5 + 0.2 * t2) - 75 * t1);

  double ZlEff, ErEff, WEff, ZlEffFreq, ErEffFreq;
  msline::analyseQuasiStatic (W1, h, t, er, SModel, ZlEff, ErEff, WEff);
  msline::analyseDispersion  (W1, h, er, ZlEff, ErEff, frequency, DModel,
			      ZlEffFreq, ErEffFreq);
  double L1 = ZlEffFreq * std::sqrt (ErEffFreq) / C0;

  msline::analyseQuasiStatic (W2, h, t, er, SModel, ZlEff, ErEff, WEff);
  msline::analyseDispersion  (W2, h, er, ZlEff, ErEff, frequency, DModel,
			      ZlEffFreq, ErEffFreq);
  double L2 = ZlEffFreq * std::sqrt (ErEffFreq) / C0;

  Ls /= (L1 + L2);
  L1 *= Ls;
  L2 *= Ls;

  // build Z-parameter matrix
  nr_complex_t z21 = nr_complex_t (0.0, -0.5e12 / (pi * frequency * Cs));
  nr_complex_t z11 = nr_complex_t (0.0, 2e-9 * pi * frequency * L1) + z21;
  nr_complex_t z22 = nr_complex_t (0.0, 2e-9 * pi * frequency * L2) + z21;
  matrix z (2);
  z.set (0, 0, z11);
  z.set (0, 1, z21);
  z.set (1, 0, z21);
  z.set (1, 1, z22);
  return z;
}

void msstep::initDC (void) {
  // a DC short (voltage source V = 0 volts)
  setVoltageSources (1);
  setInternalVoltageSource (1);
  allocMatrixMNA ();
  clearY ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void msstep::initAC (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
}

void msstep::calcAC (double frequency) {
  setMatrixY (ztoy (calcMatrixZ (frequency)));
}

void msstep::initTR (void) {
  initDC ();
}

// properties
PROP_REQ [] = {
  { "W1", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "W2", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  { "MSDispModel", PROP_STR, { PROP_NO_VAL, "Kirschning" }, PROP_RNG_DIS },
  { "MSModel", PROP_STR, { PROP_NO_VAL, "Hammerstad" }, PROP_RNG_MOD },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t msstep::cirdef =
  { "MSTEP", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
