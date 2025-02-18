/*
 * nand.cpp - logical nand class implementation
 *
 * Copyright (C) 2005, 2009 Stefan Jahn <stefan@lkcc.org>
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
#include "digital.h"
#include "nand.h"

using namespace qucs;

logicnand::logicnand () : digital () {
  type = CIR_NAND;
  setVariableSized (true);
}

void logicnand::calcOutput (void) {
  double v = getPropertyDouble ("V");
  double n = getSize () - 1;
  double x;
  for (x = 0, i = 0; i < n; i++) {
    x += 2 / (1 + calcTransfer (i));
  }
  Vout = v * (1 - n / x);
}

void logicnand::calcDerivatives (void) {
  double n = getSize () - 1;
  double x;
  for (int k = 0; k < n; k++) {
    for (x = 0, i = 0; i < n; i++) {
      x += 2 / (1 + calcTransfer (i));
    }
    x *= (1 + calcTransfer (k));
    g[k] = -2 * n * calcDerivative (k) / x / x;
  }
}

// properties
PROP_REQ [] = {
  { "V", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE }, PROP_NO_PROP };
PROP_OPT [] = {
  { "t", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "TR", PROP_REAL, { 10, PROP_NO_STR }, PROP_RNGII (1, 100) },
  PROP_NO_PROP };
struct define_t logicnand::cirdef =
  { "NAND",
    PROP_NODES, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
