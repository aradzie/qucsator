/*
 * inductor.cpp - inductor class implementation
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

#include "inductor.h"

using namespace qucs;

inductor::inductor() : circuit(2) {
  type = CIR_INDUCTOR;
  setISource(true);
}

void inductor::calcSP(double frequency) {
  double l = getPropertyDouble("L") / z0;
  nr_complex_t z = nr_complex_t(0, 2.0 * pi * frequency * l);
  setS(NODE_1, NODE_1, z / (z + 2.0));
  setS(NODE_2, NODE_2, z / (z + 2.0));
  setS(NODE_1, NODE_2, 2.0 / (z + 2.0));
  setS(NODE_2, NODE_1, 2.0 / (z + 2.0));
}

void inductor::initDC() {
  setVoltageSources(1);
  allocMatrixMNA();
  voltageSource(VSRC_1, NODE_1, NODE_2);
}

void inductor::calcDC() { clearY(); }

void inductor::initAC() {
  double l = getPropertyDouble("L");
  if (l != 0.0) {
    // for non-zero inductance usual MNA entries
    setVoltageSources(0);
    allocMatrixMNA();
  } else {
    // for zero inductance create a zero voltage source
    initDC();
    calcDC();
  }
}

void inductor::calcAC(double frequency) {
  double l = getPropertyDouble("L");
  if (l != 0.0) {
    // for non-zero inductance usual MNA entries
    nr_complex_t y = nr_complex_t(0, -1 / (2.0 * pi * frequency * l));
    setY(NODE_1, NODE_1, +y);
    setY(NODE_2, NODE_2, +y);
    setY(NODE_1, NODE_2, -y);
    setY(NODE_2, NODE_1, -y);
  }
}

void inductor::initTR() {
  initDC();
  clearY();
  setStates(2);
}

#define fState 0 // flux state
#define vState 1 // voltage state

void inductor::calcTR(double) {
  double l = getPropertyDouble("L");
  double i = real(getJ(VSRC_1));
  /* apply initial condition if requested */
  if (getMode() == MODE_INIT && isPropertyGiven("I")) {
    i = getPropertyDouble("I");
  }
  double r, v;
  setState(fState, i * l);
  integrate(fState, l, r, v);
  setD(VSRC_1, VSRC_1, -r);
  setE(VSRC_1, v);
}

void inductor::initHB() {
  setVoltageSources(1);
  setInternalVoltageSource(1);
  allocMatrixMNA();
  voltageSource(VSRC_1, NODE_1, NODE_2);
}

void inductor::calcHB(double frequency) {
  double l = getPropertyDouble("L");
  setD(VSRC_1, VSRC_1, -l * 2 * pi * frequency);
}

PROP_REQ[] = {
    {"L", PROP_REAL, {1e-9, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
PROP_OPT[] = {
    {"I", PROP_REAL, {0, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
struct define_t inductor::cirdef = {
    "L", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF,
};
