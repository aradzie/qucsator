/*
 * mscross.cpp - microstrip cross-junction class implementation
 *
 * Copyright (C) 2004, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "mscross.h"

using namespace qucs;

mscross::mscross () : circuit (6) {
  type = CIR_MSCROSS;
}

void mscross::initModel (void) {
  setNode (NODE_5, createInternal (getName (), "i13"));
  setNode (NODE_6, createInternal (getName (), "i24"));
}

void mscross::initSP (void) {
  initModel ();
  allocMatrixS ();
}

void mscross::calcSP (double frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

void mscross::initDC (void) {
  initModel ();
  setVoltageSources (5);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_5);
  voltageSource (VSRC_2, NODE_3, NODE_5);
  voltageSource (VSRC_3, NODE_2, NODE_6);
  voltageSource (VSRC_4, NODE_4, NODE_6);
  voltageSource (VSRC_5, NODE_5, NODE_6);
}

void mscross::initAC (void) {
  initModel ();
  setVoltageSources (0);
  allocMatrixMNA ();
}

void mscross::calcAC (double frequency) {
  setMatrixY (calcMatrixY (frequency));
}

double mscross::capCorrection (double W, double f) {
  substrate * subst = getSubstrate ();
  double er = subst->getPropertyDouble ("er");
  double h  = subst->getPropertyDouble ("h");
  double t  = subst->getPropertyDouble ("t");
  const char * SModel = getPropertyString ("MSModel");
  const char * DModel = getPropertyString ("MSDispModel");
  double Zl1, Er1, Zl2, Er2;
  double ZlEff, ErEff, WEff;
  msline::analyseQuasiStatic (W, h, t, 9.9, SModel, ZlEff, ErEff, WEff);
  msline::analyseDispersion  (W, h, 9.9, ZlEff, ErEff, f, DModel,
                              Zl1, Er1);
  msline::analyseQuasiStatic (W, h, t, er, SModel, ZlEff, ErEff, WEff);
  msline::analyseDispersion  (W, h, er, ZlEff, ErEff, f, DModel,
                              Zl2, Er2);
  return Zl1 / Zl2 * qucs::sqrt (Er2 / Er1);
}

double mscross::calcCap (double W1, double h, double W2) {
  double W1h = W1 / h;
  double W2h = W2 / h;
  double X = qucs::log10 (W1h) * (86.6 * W2h - 30.9 * qucs::sqrt (W2h) + 367) +
    cubic (W2h) + 74 * W2h + 130;
  return 1e-12 * W1 * (0.25 * X * qucs::pow (W1h, -1.0 / 3.0) - 60 +
			      1 / W2h / 2 - 0.375 * W1h * (1 - W2h));
 }

double mscross::calcInd (double W1, double h, double W2) {
  double W1h = W1 / h;
  double W2h = W2 / h;
  double Y = 165.6 * W2h + 31.2 * qucs::sqrt (W2h) - 11.8 * sqr (W2h);
  return 1e-9 * h * (Y * W1h - 32 * W2h + 3) * qucs::pow (W1h, -1.5);
}

matrix mscross::calcMatrixY (double f) {
  double W1 = getPropertyDouble ("W1");
  double W2 = getPropertyDouble ("W2");
  double W3 = getPropertyDouble ("W3");
  double W4 = getPropertyDouble ("W4");
  substrate * subst = getSubstrate ();
  double h  = subst->getPropertyDouble ("h");
  double W1h = (W1 + W3) / 2 / h;
  double W2h = (W2 + W4) / 2 / h;
  double C1, C2, C3, C4, L1, L2, L3, L4, L5;

  // apply asymmetric modifications of original model
  C1 = calcCap (W1, h, (W2 + W4) / 2);
  C2 = calcCap (W2, h, (W1 + W3) / 2);
  C3 = calcCap (W3, h, (W4 + W2) / 2);
  C4 = calcCap (W4, h, (W3 + W1) / 2);

  L1 = calcInd (W1, h, (W2 + W4) / 2);
  L2 = calcInd (W2, h, (W1 + W3) / 2);
  L3 = calcInd (W3, h, (W4 + W2) / 2);
  L4 = calcInd (W4, h, (W3 + W1) / 2);

  L5 = 1e-9 * h * (5 * W2h * qucs::cos (pi / 2 * (1.5 - W1h)) -
		   (1 + 7 / W1h ) / W2h - 337.5);

  // center inductance correction
  L5 = L5 * 0.8;

  // capacitance corrections
  C1 = C1 * capCorrection (W1, f);
  C2 = C2 * capCorrection (W2, f);
  C3 = C3 * capCorrection (W3, f);
  C4 = C4 * capCorrection (W4, f);

  // compute admittance matrix
  double o = 2 * pi * f;
  nr_complex_t yc1 = nr_complex_t (0, o * C1);
  nr_complex_t yc2 = nr_complex_t (0, o * C2);
  nr_complex_t yc3 = nr_complex_t (0, o * C3);
  nr_complex_t yc4 = nr_complex_t (0, o * C4);
  nr_complex_t yl1 = 1.0 / nr_complex_t (0, o * L1);
  nr_complex_t yl2 = 1.0 / nr_complex_t (0, o * L2);
  nr_complex_t yl3 = 1.0 / nr_complex_t (0, o * L3);
  nr_complex_t yl4 = 1.0 / nr_complex_t (0, o * L4);
  nr_complex_t yl5 = 1.0 / nr_complex_t (0, o * L5);
  matrix Y (6);
  Y.set (NODE_1, NODE_1, yl1 + yc1);
  Y.set (NODE_2, NODE_2, yl2 + yc2);
  Y.set (NODE_3, NODE_3, yl3 + yc3);
  Y.set (NODE_4, NODE_4, yl4 + yc4);
  Y.set (NODE_1, NODE_5, -yl1); Y.set (NODE_5, NODE_1, -yl1);
  Y.set (NODE_3, NODE_5, -yl3); Y.set (NODE_5, NODE_3, -yl3);
  Y.set (NODE_2, NODE_6, -yl2); Y.set (NODE_6, NODE_2, -yl2);
  Y.set (NODE_4, NODE_6, -yl4); Y.set (NODE_6, NODE_4, -yl4);
  Y.set (NODE_5, NODE_6, -yl5); Y.set (NODE_6, NODE_5, -yl5);
  Y.set (NODE_5, NODE_5, yl1 + yl3 + yl5);
  Y.set (NODE_6, NODE_6, yl2 + yl4 + yl5);
  return Y;
}

void mscross::initTR (void) {
  initDC ();
}

// properties
PROP_REQ [] = {
  { "W1", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "W2", PROP_REAL, { 2e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "W3", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "W4", PROP_REAL, { 2e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  { "MSDispModel", PROP_STR, { PROP_NO_VAL, "Kirschning" }, PROP_RNG_DIS },
  { "MSModel", PROP_STR, { PROP_NO_VAL, "Hammerstad" }, PROP_RNG_MOD },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t mscross::cirdef =
  { "MCROSS", 4, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
