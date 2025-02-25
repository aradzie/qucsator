/*
 * capacitor.cpp - capacitor class implementation
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

#include "capacitor.h"

using namespace qucs;

capacitor::capacitor() : circuit(2) {
  type = CIR_CAPACITOR;
  setISource(true);
}

void capacitor::calcSP(double frequency) {
  double c = getPropertyDouble("C") * z0;
  nr_complex_t y = 2.0 * nr_complex_t(0, 2.0 * pi * frequency * c);
  setS(NODE_1, NODE_1, 1.0 / (1.0 + y));
  setS(NODE_2, NODE_2, 1.0 / (1.0 + y));
  setS(NODE_1, NODE_2, y / (1.0 + y));
  setS(NODE_2, NODE_1, y / (1.0 + y));
}

void capacitor::initDC() { allocMatrixMNA(); }

void capacitor::initAC() { allocMatrixMNA(); }

void capacitor::calcAC(double frequency) {
  double c = getPropertyDouble("C");
  nr_complex_t y = nr_complex_t(0, 2.0 * pi * frequency * c);
  setY(NODE_1, NODE_1, +y);
  setY(NODE_2, NODE_2, +y);
  setY(NODE_1, NODE_2, -y);
  setY(NODE_2, NODE_1, -y);
}

#define qState 0 // charge state
#define cState 1 // current state

void capacitor::initTR() {
  setStates(2);
  initDC();
}

void capacitor::calcTR(double) {
  /* if this is a controlled capacitance then do nothing here */
  if (hasProperty("Controlled")) {
    return;
  }
  double c = getPropertyDouble("C");
  double v = real(getV(NODE_1) - getV(NODE_2));
  /* apply initial condition if requested */
  if (getMode() == MODE_INIT && isPropertyGiven("V")) {
    v = getPropertyDouble("V");
  }
  double g, i;
  setState(qState, c * v);
  integrate(qState, c, g, i);
  setY(NODE_1, NODE_1, +g);
  setY(NODE_2, NODE_2, +g);
  setY(NODE_1, NODE_2, -g);
  setY(NODE_2, NODE_1, -g);
  setI(NODE_1, -i);
  setI(NODE_2, +i);
}

void capacitor::initHB() { initAC(); }

void capacitor::calcHB(double frequency) { calcAC(frequency); }

PROP_REQ[] = {
    {"C", PROP_REAL, {1e-12, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
PROP_OPT[] = {
    {"V", PROP_REAL, {0, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
struct define_t capacitor::cirdef = {
    "C", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF,
};
