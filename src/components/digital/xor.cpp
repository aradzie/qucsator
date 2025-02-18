/*
 * xor.cpp - logical xor class implementation
 *
 * Copyright (C) 2005, 2006, 2009 Stefan Jahn <stefan@lkcc.org>
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
#include "xor.h"

using namespace qucs;

logicxor::logicxor () : digital () {
  type = CIR_XOR;
  setVariableSized (true);
}

void logicxor::calcOutput (void) {
  double v = getPropertyDouble ("V");
  double n = getSize () - 1;
  double x;
  for (x = 1, i = 0; i < n; i++) {
    x *= -calcTransferX (i);
  }
  Vout = v / 2 * (1 - x);
}

void logicxor::calcDerivatives (void) {
  double n = getSize () - 1;
  double x;
  for (int k = 0; k < n; k++) {
    for (x = 1, i = 0; i < n; i++) {
      if (i != k) x *= -calcTransferX (i);
    }
    g[k] = 0.5 * calcDerivativeX (k) * x;
  }
}

// properties
PROP_REQ [] = {
  { "V", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE }, PROP_NO_PROP };
PROP_OPT [] = {
  { "t", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "TR", PROP_REAL, { 10, PROP_NO_STR }, PROP_RNGII (1, 100) },
  PROP_NO_PROP };
struct define_t logicxor::cirdef =
  { "XOR",
    PROP_NODES, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
