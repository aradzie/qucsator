/*
 * inoise.cpp - noise current source class implementation
 *
 * Copyright (C) 2004, 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "inoise.h"

using namespace qucs;

inoise::inoise () : circuit (2) {
  type = CIR_INOISE;
}

void inoise::initSP (void) {
  allocMatrixS ();
  setS (NODE_1, NODE_1, 1.0);
  setS (NODE_1, NODE_2, 0.0);
  setS (NODE_2, NODE_1, 0.0);
  setS (NODE_2, NODE_2, 1.0);
}

void inoise::calcNoiseSP (double frequency) {
  double i = getPropertyDouble ("i");
  double e = getPropertyDouble ("e");
  double c = getPropertyDouble ("c");
  double a = getPropertyDouble ("a");
  double ipsd = i / (a + c * qucs::pow (frequency, e)) / kB / T0 * z0;
  setN (NODE_1, NODE_1, +ipsd); setN (NODE_2, NODE_2, +ipsd);
  setN (NODE_1, NODE_2, -ipsd); setN (NODE_2, NODE_1, -ipsd);
}

void inoise::calcNoiseAC (double frequency) {
  double i = getPropertyDouble ("i");
  double e = getPropertyDouble ("e");
  double c = getPropertyDouble ("c");
  double a = getPropertyDouble ("a");
  double ipsd = i / (a + c * qucs::pow (frequency, e)) / kB / T0;
  setN (NODE_1, NODE_1, +ipsd); setN (NODE_2, NODE_2, +ipsd);
  setN (NODE_1, NODE_2, -ipsd); setN (NODE_2, NODE_1, -ipsd);
}

// properties
PROP_REQ [] = {
  { "i", PROP_REAL, { 1e-6, PROP_NO_STR }, PROP_POS_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "a", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "c", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "e", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  PROP_NO_PROP };
struct define_t inoise::cirdef =
  { "Inoise", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
