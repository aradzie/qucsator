/*
 * resistor.cpp - resistor class implementation
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include "resistor.h"

using namespace qucs;

resistor::resistor() : circuit(2) { type = CIR_RESISTOR; }

void resistor::initSP() {
  initModel();
  allocMatrixS();
}

void resistor::calcSP(double) {
  // calculate S-parameters
  const double z = getScaledProperty("R") / z0;
  setS(NODE_1, NODE_1, z / (z + 2));
  setS(NODE_2, NODE_2, z / (z + 2));
  setS(NODE_1, NODE_2, 2 / (z + 2));
  setS(NODE_2, NODE_1, 2 / (z + 2));
}

void resistor::calcNoiseSP(double) {
  // calculate noise correlation matrix
  const double r = getScaledProperty("R");
  const double T = getPropertyDouble("Temp");
  const double f = celsius2kelvin(T) * 4.0 * r * z0 / sqr(2.0 * z0 + r) / T0;
  setN(NODE_1, NODE_1, +f);
  setN(NODE_2, NODE_2, +f);
  setN(NODE_1, NODE_2, -f);
  setN(NODE_2, NODE_1, -f);
}

void resistor::calcNoiseAC(double) {
  // calculate noise current correlation matrix
  const double r = getScaledProperty("R");
  if (r > 0.0 || r < 0.0) {
    const double T = getPropertyDouble("Temp");
    const double f = celsius2kelvin(T) / T0 * 4.0 / r;
    setN(NODE_1, NODE_1, +f);
    setN(NODE_2, NODE_2, +f);
    setN(NODE_1, NODE_2, -f);
    setN(NODE_2, NODE_1, -f);
  }
}

void resistor::initModel() {
  /* if this is a controlled resistor then do nothing here */
  if (hasProperty("Controlled"))
    return;

  double T = getPropertyDouble("Temp");
  double Tn = getPropertyDouble("Tnom");
  double R = getPropertyDouble("R");
  double DT = T - Tn;

  // compute R temperature dependency
  double Tc1 = getPropertyDouble("Tc1");
  double Tc2 = getPropertyDouble("Tc2");
  R = R * (1 + DT * (Tc1 + Tc2 * DT));
  setScaledProperty("R", R);
}

void resistor::initDC() {
  initModel();

  const double r = getScaledProperty("R");
  if (r != 0.0) {
    // for non-zero resistances usual MNA entries
    setVoltageSources(0);
    allocMatrixMNA();
    const double g = 1.0 / r;
    setY(NODE_1, NODE_1, +g);
    setY(NODE_2, NODE_2, +g);
    setY(NODE_1, NODE_2, -g);
    setY(NODE_2, NODE_1, -g);
  } else {
    // for zero resistances create a zero voltage source
    setVoltageSources(1);
    setInternalVoltageSource(1);
    allocMatrixMNA();
    voltageSource(VSRC_1, NODE_1, NODE_2);
  }
}

/* The calcDC() function is here partly implemented again because the
   circuit can be used to simulate controlled non-zero resistances. */
void resistor::calcDC() {
  const double r = getScaledProperty("R");
  if (r != 0.0) {
    // for non-zero resistances usual MNA entries
    const double g = 1.0 / r;
    setY(NODE_1, NODE_1, +g);
    setY(NODE_2, NODE_2, +g);
    setY(NODE_1, NODE_2, -g);
    setY(NODE_2, NODE_1, -g);
  }
}

void resistor::initAC() { initDC(); }

void resistor::calcAC(double) { calcDC(); }

void resistor::initTR() { initDC(); }

void resistor::calcTR(double) { calcDC(); }

// Initialize computation of MNA matrix entries for HB.
void resistor::initHB() {
  initModel();
  double r = getScaledProperty("R");
  setVoltageSources(1);
  setInternalVoltageSource(1);
  allocMatrixMNA();
  voltageSource(VSRC_1, NODE_1, NODE_2);
  setD(VSRC_1, VSRC_1, -r);
}

// properties
PROP_REQ[] = {
    {"R", PROP_REAL, {50, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
PROP_OPT[] = {
    {"Temp", PROP_REAL, {26.85, PROP_NO_STR}, PROP_MIN_VAL(K)},
    {"Tc1", PROP_REAL, {0, PROP_NO_STR}, PROP_NO_RANGE},
    {"Tc2", PROP_REAL, {0, PROP_NO_STR}, PROP_NO_RANGE},
    {"Tnom", PROP_REAL, {26.85, PROP_NO_STR}, PROP_MIN_VAL(K)},
    PROP_NO_PROP,
};
struct define_t resistor::cirdef = {
    "R", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF,
};
