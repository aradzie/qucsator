/*
 * trsolver.cpp - transient solver class implementation
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2009 Stefan Jahn <stefan@lkcc.org>
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

#include "circuit.h"
#include "exception.h"
#include "exceptionstack.h"
#include "logging.h"
#include "trsolver.h"

#define STEPDEBUG 0   // set to zero for release
#define BREAKPOINTS 0 // exact breakpoint calculation

namespace qucs {


// Performs the initial DC analysis.
int trsolver::dcAnalysis() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::dcAnalysis()\n", getName());

  int error = 0;

  // First calculate the initial state using the non-linear DC analysis.
  setDescription("initial DC");
  initDC();
  setCalculation((calculate_func_t)&calcDC);
  solve_pre();
  applyNodeset();

  // Run the DC solver once.
  error = solve_nonlinear();

  if (estack.top()) {
    switch (estack.top()->getCode()) {
    case EXCEPTION_NO_CONVERGENCE:
      estack.pop();
      convHelper = CONV_LineSearch;
      applyNodeset();
      error = solve_nonlinear();
      break;
    default:
      // Otherwise return.
      estack.print();
      return -1;
    }
  }

  // Save the DC solution.
  storeDcSolution();

  // Cleanup nodal analysis solver.
  solve_post();

  if (error) {
    logprint(LOG_ERROR, "ERROR: %s: %s analysis failed\n", getName(), getDescription().c_str());
  }

  return error;
}

/* Goes through the list of circuit objects and runs its initDC() function. */
void trsolver::initDC() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::initDC()\n", getName());

  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->initDC();
  }
}

/* Goes through the list of circuit objects and runs its calcDC() function. */
void trsolver::calcDC(trsolver *self) {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::calcDC()\n", self->getName());

  circuit *root = self->getNet()->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->calcDC();
  }
}

// Stores the DC solution (node voltages and branch currents).
void trsolver::storeDcSolution() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::storeDcSolution()\n", getName());

  // cleanup solution previously
  dcSolution.clear();
  const int N = countNodes();
  const int M = countVoltageSources();
  // store all nodes except reference node
  for (int r = 0; r < N; r++) {
    struct nodelist_t *n = nlist->getNode(r);
    double gr = x->get(r);
    naentry entry(gr, 0);
    dcSolution.insert({{n->name, entry}});
    logprint(LOG_STATUS, "NOTIFY: %s: save solution entry %s=%e\n", getName(), n->name.c_str(), gr);
  }
  // store all branch currents of voltage sources
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    int vn = r - vs->getVoltageSource() + 1;
    double xg = x->get(r + N);
    naentry entry(xg, vn);
    dcSolution.insert({{vs->getName(), entry}});
    logprint(LOG_STATUS, "NOTIFY: %s: save solution entry %s=%e\n", getName(), vs->getName(), xg);
  }
}

// Recalls the DC solution (node voltages and branch currents).
void trsolver::recallDcSolution() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::recallDcSolution()\n", getName());

  const int N = countNodes();
  const int M = countVoltageSources();
  // store all nodes except reference node
  for (int r = 0; r < N; r++) {
    struct nodelist_t *n = nlist->getNode(r);
    auto na = dcSolution.find(n->name);
    if (na != dcSolution.end())
      if ((*na).second.current == 0)
        x->set(r, (*na).second.value);
  }
  // store all branch currents of voltage sources
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    int vn = r - vs->getVoltageSource() + 1;
    auto na = dcSolution.find(vs->getName());
    if (na != dcSolution.end())
      if ((*na).second.current == vn)
        x->set(r + N, (*na).second.value);
  }
}

} // namespace qucs
