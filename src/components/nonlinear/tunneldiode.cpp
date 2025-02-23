/*
 * tunneldiode.cpp - resonance tunnel diode class implementation
 *
 * Copyright (C) 2011 Michael Margraf <michael.margraf@alumni.tu-berlin.de>
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
#include "tunneldiode.h"

#define NODE_A1 0 /* cathode node */
#define NODE_A2 1 /* anode node   */

using namespace qucs;
using namespace qucs::device;

// Constructor for the diode.
tunneldiode::tunneldiode () : circuit (2) {
  type = CIR_TUNNELDIODE;
}

// Callback for initializing the DC analysis.
void tunneldiode::initDC (void) {
  // allocate MNA matrices
  allocMatrixMNA ();
}

// Calculate one branch of the tunnel current.
void tunneldiode::calcId (double U, double& I, double& G) {
  double eta  = getPropertyDouble ("eta");
  double Wr   = getPropertyDouble ("Wr");
  double dv   = getPropertyDouble ("dv");
  double de   = getPropertyDouble ("de");
  double dW   = getPropertyDouble ("dW");

  U   = Wr - Q_e*U/dv;
  de *= kB * celsius2kelvin (getPropertyDouble ("Temp"));

  double a = pi_over_2 + qucs::atan ( U / dW );

  double e = (eta - U) / de;
  double b = e;
  if (e < 15.0)  // avoid numerical overflow
    b = qucs::log (1.0 + qucs::exp ( e ));

  // current
  I = b * a;

  // derivative
  G = Q_e / dv / de / (1.0 + qucs::exp(-e)) * a - b * Q_e / dv / dW / (1.0 + sqr (U/dW));
}

// Callback for DC analysis.
void tunneldiode::calcDC (void) {
  // get device properties
  double Ip   = getPropertyDouble ("Ip");
  double A    = getPropertyDouble ("Area");
  double Tmax = getPropertyDouble ("Tmax");
  double de   = getPropertyDouble ("de");
  double eta  = getPropertyDouble ("eta");
  double Iv   = getPropertyDouble ("Iv");
  double Vv   = getPropertyDouble ("Vv");
  double nv   = getPropertyDouble ("nv");
  double T    = kB * celsius2kelvin (getPropertyDouble ("Temp"));

  // diode voltage
  Ud = real (getV (NODE_A1) - getV (NODE_A2));

  // bi-directional tunnel current
  double Ipos, Ineg, Gpos, Gneg;
  gd = Id = A * Ip * Tmax * de * T / eta / pi_over_2;
  calcId ( Ud, Ipos, Gpos);
  calcId (-Ud, Ineg, Gneg);
  Id *= Ipos - Ineg;
  gd *= Gpos + Gneg;

  // thermal-ionic current
  nv *= T / Q_e;
  double c = A * Iv / qucs::sinh (Vv / nv);
  Id += c * qucs::sinh (Ud / nv);
  gd += c * qucs::cosh (Ud / nv) / nv;

  double Ieq = Id - Ud * gd;

  // fill in I-Vector
  setI (NODE_A2, +Ieq);
  setI (NODE_A1, -Ieq);

  // fill in G-Matrix
  setY (NODE_A1, NODE_A1, +gd); setY (NODE_A2, NODE_A2, +gd);
  setY (NODE_A1, NODE_A2, -gd); setY (NODE_A2, NODE_A1, -gd);
}

// Saves operating points (voltages).
void tunneldiode::saveOperatingPoints (void) {
  double Vd = real (getV (NODE_A1) - getV (NODE_A2));
  setOperatingPoint ("Vd", Vd);
}

// Loads operating points (voltages).
void tunneldiode::loadOperatingPoints (void) {
  Ud = getOperatingPoint ("Vd");
}

// Calculates and saves operating points.
void tunneldiode::calcOperatingPoints (void) {
  double A   = getPropertyDouble ("Area");
  double Cj0 = getPropertyDouble ("Cj0");
  double M   = getScaledProperty ("M");
  double Vj  = getScaledProperty ("Vj");
  double te  = getScaledProperty ("te");

  // calculate capacitances and charges
  double Cd;

  // depletion capacitance
  double c = 1.0 + fabs(Ud) / Vj;
  Cd = A * Cj0 / qucs::pow (c, M);
  Qd = A * Cj0 * Vj / (1.0-M) * (1.0 - qucs::pow (c, 1.0 - M));

  // quantum well (diffusion) capacitance (negative because of NDR region)
  Cd -= te * gd;
  Qd -= te * Id;

  // save operating points
  setOperatingPoint ("gd", gd);
  setOperatingPoint ("Id", Id);
  setOperatingPoint ("Cd", Cd);
}

// Callback for initializing the AC analysis.
void tunneldiode::initAC (void) {
  initDC ();
}

// Build admittance matrix for AC and SP analysis.
matrix tunneldiode::calcMatrixY (double frequency) {
  double gd = getOperatingPoint ("gd");
  double Cd = getOperatingPoint ("Cd");
  nr_complex_t yd = nr_complex_t (gd, Cd * 2.0 * pi * frequency);
  matrix y (2);
  y.set (NODE_A1, NODE_A1, +yd);
  y.set (NODE_A2, NODE_A2, +yd);
  y.set (NODE_A1, NODE_A2, -yd);
  y.set (NODE_A2, NODE_A1, -yd);
  return y;
}

// Callback for the AC analysis.
void tunneldiode::calcAC (double frequency) {
  setMatrixY (calcMatrixY (frequency));
}

// Callback for S-parameter analysis.
void tunneldiode::calcSP (double frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

#define qState 0 // charge state
#define cState 1 // current state

// Callback for initializing the TR analysis.
void tunneldiode::initTR (void) {
  setStates (2);
  initDC ();
}

// Callback for the TR analysis.
void tunneldiode::calcTR (double) {
  calcDC ();

  saveOperatingPoints ();
  loadOperatingPoints ();
  calcOperatingPoints ();

  double Cd = getOperatingPoint ("Cd");
  transientCapacitance (qState, NODE_A1, NODE_A2, Cd, Ud, Qd);
}

// properties
PROP_REQ [] = {
  { "Ip", PROP_REAL, { 4.0e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Iv", PROP_REAL, { 0.6e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vv", PROP_REAL, { 0.8, PROP_NO_STR }, PROP_POS_RANGE },

  { "Cj0", PROP_REAL, { 80e-15, PROP_NO_STR }, PROP_POS_RANGE },
  { "M", PROP_REAL, { 0.5, PROP_NO_STR }, PROP_RNGII (0, 2) },
  { "Vj", PROP_REAL, { 0.5, PROP_NO_STR }, PROP_RNGXI (0, 10) },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Wr", PROP_REAL, { 2.7e-20, PROP_NO_STR }, PROP_POS_RANGE },
  { "eta", PROP_REAL, { 1e-20, PROP_NO_STR }, PROP_POS_RANGE },
  { "dW", PROP_REAL, { 4.5e-21, PROP_NO_STR }, PROP_POS_RANGE },
  { "Tmax", PROP_REAL, { 0.95, PROP_NO_STR }, PROP_POS_RANGE },
  { "de", PROP_REAL, { 0.9, PROP_NO_STR }, PROP_POS_RANGE },
  { "dv", PROP_REAL, { 2.0, PROP_NO_STR }, PROP_POS_RANGE },
  { "nv", PROP_REAL, { 16, PROP_NO_STR }, PROP_POS_RANGE },
  { "te", PROP_REAL, { 0.6e-12, PROP_NO_STR }, PROP_POS_RANGE },
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  { "Area", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGEX },
  PROP_NO_PROP };
struct define_t tunneldiode::cirdef =
  { "RTD", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
