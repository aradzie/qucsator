/*
 * diac.cpp - diac class implementation
 *
 * Copyright (C) 2008 Stefan Jahn <stefan@lkcc.org>
 * Copyright (C) 2008 Michael Margraf <Michael.Margraf@alumni.TU-Berlin.DE>
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
#include "device.h"
#include "devstates.h"
#include "diac.h"

#define NODE_A1 0 /* anode 1 */
#define NODE_A2 1 /* anode 2 (cathode) */
#define NODE_IN 2 /* internal node */

using namespace qucs;
using namespace qucs::device;

// Constructor for the diac.
diac::diac () : circuit (3) {
  type = CIR_DIAC;
}

// Callback for initializing the DC analysis.
void diac::initDC (void) {
  Ud_last = 0.0;
  // allocate MNA matrices
  allocMatrixMNA ();
  // create internal node
  setInternalNode (NODE_IN, "int");
}

// Callback for DC analysis.
void diac::calcDC (void) {
  calcTheModel (false);
}

void diac::calcTheModel (bool last) {
  // get device properties
  double Ubo = getPropertyDouble ("Vbo");
  double Ibo = getPropertyDouble ("Ibo");
  double Is  = getPropertyDouble ("Is");
  double N   = getPropertyDouble ("N");
  double gi  = 1.0 / getPropertyDouble ("Ri");
  double T   = getPropertyDouble ("Temp");

  bool isOn;
  if (last)
    Ud = fabs (Ud_last);
  else
    Ud = fabs (real (getV (NODE_A1) - getV (NODE_IN)));
  isOn = Ud > (Ibo / gi);

  double Ut, Ieq, Vd;

  if (isOn)
    Ut = N * celsius2kelvin (T) * kBoverQ;
  else
    Ut  = Ubo / std::log (Ibo / Is);

  Vd = Ud = real (getV (NODE_IN) - getV (NODE_A2));
  Ud = fabs (Ud) / Ut;
  Id = sign (Vd) * Is;

  if (Ud >= 80.0) {
    Id *= std::exp (80.0) * (1.0 + Ud - 80.0) - 1.0;
    Ud  = 80.0;
  }
  else
    Id *= std::exp (Ud) - 1.0;

  gd  = Is / Ut * std::exp (Ud);
  Ieq = Id - Vd * gd;

  // fill in I-Vector
  setI (NODE_A2, +Ieq);
  setI (NODE_IN, -Ieq);
  setI (NODE_A1, +0.0);

  // fill in G-Matrix
  setY (NODE_A2, NODE_A2, +gd); setY (NODE_IN, NODE_IN, +gd);
  setY (NODE_A2, NODE_IN, -gd); setY (NODE_IN, NODE_A2, -gd);
  setY (NODE_A1, NODE_A1, +gi); addY (NODE_IN, NODE_IN, +gi);
  setY (NODE_A1, NODE_IN, -gi); setY (NODE_IN, NODE_A1, -gi);
}

// Saves operating points (voltages).
void diac::saveOperatingPoints (void) {
  double Vd = real (getV (NODE_IN) - getV (NODE_A2));
  double Vi = real (getV (NODE_A1) - getV (NODE_IN));
  setOperatingPoint ("Vd", Vd);
  setOperatingPoint ("Vi", Vi);
}

// Loads operating points (voltages).
void diac::loadOperatingPoints (void) {
  Ud = getOperatingPoint ("Vd");
  Ui = getOperatingPoint ("Vi");
}

// Calculates and saves operating points.
void diac::calcOperatingPoints (void) {
  double Cj0 = getPropertyDouble ("Cj0");
  // calculate capacitances and charges
  double Ci;
  Ci = Cj0;
  Qi = Cj0 * Ud;
  // save operating points
  setOperatingPoint ("gi", gi);
  setOperatingPoint ("gd", gd);
  setOperatingPoint ("Id", Id);
  setOperatingPoint ("Ci", Ci);
}

// Callback for initializing the AC analysis.
void diac::initAC (void) {
  initDC ();
}

// Build admittance matrix for AC and SP analysis.
matrix diac::calcMatrixY (double frequency) {
  double gd = getOperatingPoint ("gd");
  double gi = getOperatingPoint ("gi");
  double Ci = getOperatingPoint ("Ci");
  nr_complex_t yd = nr_complex_t (gd, Ci * 2.0 * pi * frequency);
  matrix y (3);
  y.set (NODE_A2, NODE_A2, +yd);
  y.set (NODE_IN, NODE_IN, +yd +gi);
  y.set (NODE_A2, NODE_IN, -yd);
  y.set (NODE_IN, NODE_A2, -yd);
  y.set (NODE_A1, NODE_A1, +gi);
  y.set (NODE_A1, NODE_IN, -gi);
  y.set (NODE_IN, NODE_A1, -gi);
  return y;
}

// Callback for the AC analysis.
void diac::calcAC (double frequency) {
  setMatrixY (calcMatrixY (frequency));
}

// Callback for S-parameter analysis.
void diac::calcSP (double frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

#define qState 0 // charge state
#define cState 1 // current state

// Callback for initializing the TR analysis.
void diac::initTR (void) {
  setStates (2);
  initDC ();
  time_prev = -1.0;
}

// Callback for the TR analysis.
void diac::calcTR (double time) {
  if (time_prev < time) {
    time_prev = time;
    Ud_last = real (getV (NODE_A1) - getV (NODE_IN));
  }
  calcTheModel (true);

  saveOperatingPoints ();
  loadOperatingPoints ();
  calcOperatingPoints ();

  double Ci = getOperatingPoint ("Ci");
  transientCapacitance (qState, NODE_IN, NODE_A2, Ci, Ud, Qi);
}

// properties
PROP_REQ [] = {
  { "Ibo", PROP_REAL, { 50e-6, PROP_NO_STR }, PROP_POS_RANGEX },
  { "Vbo", PROP_REAL, { 30, PROP_NO_STR }, PROP_POS_RANGEX },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Cj0", PROP_REAL, { 10e-12, PROP_NO_STR }, PROP_POS_RANGE },
  { "Is", PROP_REAL, { 1e-10, PROP_NO_STR }, PROP_POS_RANGE },
  { "N", PROP_REAL, { 2.0, PROP_NO_STR }, PROP_RNGII (0.1, 100) },
  { "Ri", PROP_REAL, { 10.0, PROP_NO_STR }, PROP_POS_RANGEX },
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  PROP_NO_PROP };
struct define_t diac::cirdef =
  { "Diac", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
