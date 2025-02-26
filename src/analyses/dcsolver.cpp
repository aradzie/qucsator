/*
 * dcsolver.cpp - DC solver class implementation
 *
 * Copyright (C) 2003-2008 Stefan Jahn <stefan@lkcc.org>
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

#include <cstring>

#include "analysis.h"
#include "circuit.h"
#include "constants.h"
#include "dcsolver.h"
#include "exception.h"
#include "exceptionstack.h"
#include "logging.h"
#include "nasolver.h"
#include "net.h"
#include "netdefs.h"

namespace qucs {

constexpr int helpers[] = {
    CONV_SourceStepping,  //
    CONV_GMinStepping,    //
    CONV_SteepestDescent, //
    CONV_LineSearch,      //
    CONV_Attenuation,     //
    -1,
};

dcsolver::dcsolver() : nasolver<double>() {
  saveOPs = 0;
  type = ANALYSIS_DC;
  setDescription("DC");
}

dcsolver::dcsolver(char *n) : nasolver<double>(n) {
  saveOPs = 0;
  type = ANALYSIS_DC;
  setDescription("DC");
}

dcsolver::~dcsolver() {}

int dcsolver::solve() {
  logprint(LOG_STATUS, "NOTIFY: %s: dcsolver::solve()\n", getName());

  saveOPs |= !strcmp(getPropertyString("saveOPs"), "yes") ? SAVE_OPS : 0;
  saveOPs |= !strcmp(getPropertyString("saveAll"), "yes") ? SAVE_ALL : 0;

  // initialize node voltages, first guess for non-linear circuits and
  // generate extra circuits if necessary
  initDC();
  setCalculation((calculate_func_t)&calcDC);

  const char *const solver = getPropertyString("Solver");
  if (!strcmp(solver, "CroutLU"))
    eqnAlgo = ALGO_LU_DECOMPOSITION_CROUT;
  else if (!strcmp(solver, "DoolittleLU"))
    eqnAlgo = ALGO_LU_DECOMPOSITION_DOOLITTLE;
  else if (!strcmp(solver, "HouseholderQR"))
    eqnAlgo = ALGO_QR_DECOMPOSITION;
  else if (!strcmp(solver, "HouseholderLQ"))
    eqnAlgo = ALGO_QR_DECOMPOSITION_LS;
  else if (!strcmp(solver, "GolubSVD"))
    eqnAlgo = ALGO_SV_DECOMPOSITION;

  // Allocate the nodes, SLE matrices and vectors.
  solve_pre();

  convHelper = CONV_None;

  int error = 0;
  if (!subnet->isNonLinear()) {
    // Start the linear solver.
    error = solve_linear();
  } else {
    int helperIndex = 0;

    // Iterate to find the solution.
    bool retry;
    do {
      retry = false;
      // Use the user provided node voltages.
      applyNodeset();
      // Run the DC solver once.
      error = solve_nonlinear();
#if DEBUG
      if (!error) {
        logprint(LOG_STATUS, "NOTIFY: %s: convergence reached after %d iterations\n", getName(),
                 iterations);
      }
#endif /* DEBUG */
      if (estack.top()) {
        switch (estack.top()->getCode()) {
        case EXCEPTION_NO_CONVERGENCE:
          estack.pop();
          convHelper = helpers[helperIndex++];
          if (convHelper != -1) {
            logprint(LOG_ERROR,
                     "WARNING: %s: %s analysis failed, using fallback "
                     "#%d (%s)\n",
                     getName(), getDescription().c_str(), helperIndex, getHelperDescription());
            retry = true;
            restartDC();
          }
          break;
        default:
          estack.print();
          return -1;
        }
      }
    } while (retry);
  }

  // save results and cleanup the solver
  saveOperatingPoints();
  saveResults("V", "I", saveOPs);

  // Free the nodes, SLE matrices and vectors.
  solve_post();

  return error;
}

/* Goes through the list of circuit objects and runs their initDC() function. */
void dcsolver::initDC() {
  logprint(LOG_STATUS, "NOTIFY: %s: dcsolver::initDC()\n", getName());

  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->initDC();
  }
}

/* Goes through the list of circuit objects and runs their calcDC() function. */
void dcsolver::calcDC(dcsolver *self) {
  logprint(LOG_STATUS, "NOTIFY: %s: dcsolver::calcDC()\n", self->getName());

  circuit *root = self->getNet()->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->calcDC();
  }
}

/* Goes through the list of non-linear circuit objects
 * and runs its saveOperatingPoints() function. */
void dcsolver::saveOperatingPoints() {
  logprint(LOG_STATUS, "NOTIFY: %s: dcsolver::saveOperatingPoints()\n", getName());

  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->isNonLinear()) {
      c->saveOperatingPoints();
    }
  }
}

PROP_REQ[] = {PROP_NO_PROP};
PROP_OPT[] = {
    {"MaxIter", PROP_INT, {150, PROP_NO_STR}, PROP_RNGII(2, 10000)},
    {"abstol", PROP_REAL, {1e-12, PROP_NO_STR}, PROP_RNG_X01I},
    {"vntol", PROP_REAL, {1e-6, PROP_NO_STR}, PROP_RNG_X01I},
    {"reltol", PROP_REAL, {1e-3, PROP_NO_STR}, PROP_RNG_X01I},
    {"saveOPs", PROP_STR, {PROP_NO_VAL, "no"}, PROP_RNG_YESNO},
    {"Temp", PROP_REAL, {26.85, PROP_NO_STR}, PROP_MIN_VAL(K)},
    {"saveAll", PROP_STR, {PROP_NO_VAL, "no"}, PROP_RNG_YESNO},
    {"Solver", PROP_STR, {PROP_NO_VAL, "CroutLU"}, PROP_RNG_SOL},
    PROP_NO_PROP,
};
struct define_t dcsolver::anadef = {"DC", 0, PROP_ACTION, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF};

} // namespace qucs
