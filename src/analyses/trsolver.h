/*
 * trsolver.h - transient solver class definitions
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __TRSOLVER_H__
#define __TRSOLVER_H__

#include <string>
#include <unordered_map>

#include "nasolver.h"

namespace qucs {

class sweep;
class circuit;
class history;

class naentry {
public:
  naentry() = default;
  naentry(const naentry &) = default;
  naentry(const double &v, const int c) : value(v), current(c) {}
  ~naentry() = default;

public:
  double value; // MNA vector x entry.
  int current;  // Voltage source index in a circuit.
};

class trsolver final : public nasolver<double> {
public:
  ACREATOR(trsolver);
  explicit trsolver(const std::string &name);
  trsolver(const trsolver &) = delete;
  ~trsolver() override;
  int solve() override;

private:
  int dcAnalysis();
  void predictor();
  int corrector();
  void nextStates();
  void fillStates();
  void setMode(int);
  void adjustDelta(double);
  void adjustOrder(int reduce = 0);
  void initDC();
  static void calcDC(trsolver *);
  void initTR();
  static void calcTR(trsolver *);
  void saveAllResults(double);
  double checkDelta();
  void initHistory(double);
  void updateHistory(double);
  void saveHistory(circuit *);
  void predictBashford();
  void predictEuler();
  void predictGear();

  void storeDcSolution();
  void recallDcSolution();

private:
  sweep *swp;
  tvector<double> *solution[8]; // The list of previous solution vectors X.
  double predCoeff[8];          // This array is shared with the circuits (integrators).
  double corrCoeff[8];          // This array is shared with the circuits (integrators).
  double deltas[8];             // This array is shared with the circuits (integrators).
  double delta;
  double deltaMax;
  double deltaMin;
  double deltaOld;
  double stepDelta;
  int corrMethod0;  // user specified corrector method
  int corrMethod;   // current corrector method
  int predMethod;   // current predictor method
  int corrMaxOrder; // maximum corrector order
  int predMaxOrder; // maximum predictor order
  int corrOrder;    // current corrector order
  int predOrder;    // current predictor order
  int rejected;
  int converged;
  double current;
  int statSteps;
  int statRejected;
  int statIterations;
  int statConvergence;
  history *hist;

  std::unordered_map<
      /* node or circuit name */ std::string,
      /* MNA vector x entry */ naentry>
      dcSolution;
};

} // namespace qucs

#endif /* __TRSOLVER_H__ */
