/*
 * acsolver.cpp - AC solver class implementation
 *
 * Copyright (C) 2004, 2005, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include "acsolver.h"
#include "analysis.h"
#include "circuit.h"
#include "complex.h"
#include "constants.h"
#include "dataset.h"
#include "nasolver.h"
#include "net.h"
#include "netdefs.h"
#include "node.h"
#include "object.h"
#include "sweep.h"
#include "vector.h"

namespace qucs {

acsolver::acsolver() : nasolver<nr_complex_t>() {
  swp = nullptr;
  type = ANALYSIS_AC;
  setDescription("AC");
  xn = nullptr;
  noise = 0;
}

acsolver::acsolver(char *n) : nasolver<nr_complex_t>(n) {
  swp = nullptr;
  type = ANALYSIS_AC;
  setDescription("AC");
  xn = nullptr;
  noise = 0;
}

acsolver::~acsolver() {
  delete swp;
  delete xn;
}

int acsolver::solve() {
  runs++;

  noise = !strcmp(getPropertyString("Noise"), "yes") ? 1 : 0;

  if (swp == nullptr) {
    swp = createSweep("acfrequency");
  }

  // initialize node voltages, first guess for non-linear circuits and
  // generate extra circuits if necessary
  initAC();
  setCalculation((calculate_func_t)&calcAC);

  solve_pre();

  swp->reset();
  for (int i = 0; i < swp->getSize(); i++) {
    freq = swp->next();

#if DEBUG
    logprint(LOG_STATUS, "NOTIFY: %s: solving netlist for f = %e\n", getName(), freq);
#endif

    eqnAlgo = ALGO_LU_DECOMPOSITION;
    solve_linear();
    if (noise) {
      solve_noise();
    }

    saveAllResults(freq);
  }

  solve_post();

  return 0;
}

/* Goes through the list of circuit objects and runs its initAC() function. */
void acsolver::initAC() {
  logprint(LOG_STATUS, "NOTIFY: %s: acsolver::initAC()\n", getName());

  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->isNonLinear()) {
      c->calcOperatingPoints();
    }
    c->initAC();
    if (noise) {
      c->initNoiseAC();
    }
  }
}

/* Goes through the list of circuit objects and runs its calcAC() function. */
void acsolver::calcAC(acsolver *self) {
  logprint(LOG_STATUS, "NOTIFY: %s: acsolver::calcAC()\n", self->getName());

  circuit *root = self->getNet()->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->calcAC(self->freq);
    if (self->noise) {
      c->calcNoiseAC(self->freq);
    }
  }
}

/* Saves the results of a single solve() functionality (for the given frequency)
 * into the output dataset. */
void acsolver::saveAllResults(double freq) {
  logprint(LOG_STATUS, "NOTIFY: %s: acsolver::saveAllResults(%e)\n", getName(), freq);

  qucs::vector *f = data->findDependency("acfrequency");
  // add current frequency to the dependency of the output dataset
  if (f == nullptr) {
    data->addDependency(f = new qucs::vector("acfrequency"));
  }
  if (runs == 1) {
    f->add(freq);
  }
  saveResults("v", "i", 0, f);

  // additionally save noise results if requested
  if (noise) {
    saveNoiseResults(f);
  }
}

/* Computes the final noise results and puts them into the output dataset. */
void acsolver::saveNoiseResults(qucs::vector *f) {
  const int N = countNodes();
  const int M = countVoltageSources();
  for (int r = 0; r < N + M; r++) {
    // renormalise the results
    x->set(r, fabs(xn->get(r) * sqrt(kB * T0)));
  }

  // apply probe data
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (!c->isProbe())
      continue;
    int np, nn;
    double vp, vn;
    np = getNodeNr(c->getNode(NODE_1)->getName());
    vp = np > 0 ? xn->get(np - 1) : 0.0;
    nn = getNodeNr(c->getNode(NODE_2)->getName());
    vn = nn > 0 ? xn->get(nn - 1) : 0.0;
    c->setOperatingPoint("Vr", fabs((vp - vn) * sqrt(kB * T0)));
    c->setOperatingPoint("Vi", 0.0);
  }

  saveResults("vn", "in", 0, f);
}

/* Tuns the AC noise analysis.  It saves its results in
   the 'xn' vector. */
void acsolver::solve_noise() {
  const int N = countNodes();
  const int M = countVoltageSources();

  // save usual AC results
  tvector<nr_complex_t> xsave = *x;

  // create the Cy matrix
  createNoiseMatrix();

  // create noise result vector if necessary
  if (xn == nullptr) {
    xn = new tvector<double>(N + M);
  }

  // temporary result vector for transimpedances
  tvector<nr_complex_t> zn = tvector<nr_complex_t>(N + M);

  // create the MNA matrix once again and LU decompose the adjoint matrix
  createMatrix();
  A->transpose();
  eqnAlgo = ALGO_LU_FACTORIZATION_CROUT;
  solveLinearEquations();

  // ensure skipping LU decomposition
  updateMatrix = 0;
  convHelper = CONV_None;
  eqnAlgo = ALGO_LU_SUBSTITUTION_CROUT;

  // compute noise voltage for each node (and voltage source)
  for (int i = 0; i < N + M; i++) {
    z->set(0);
    z->set(i, -1); // modify right hand side appropriately
    solveLinearEquations();      // solve
    zn = *x;       // save transimpedance vector

    // compute actual noise voltage
    xn->set(i, sqrt(real(scalar(zn * (*C), conj(zn)))));
  }

  // restore usual AC results
  *x = xsave;
}

PROP_REQ[] = {
    {"Type", PROP_STR, {PROP_NO_VAL, "lin"}, PROP_RNG_TYP},
    PROP_NO_PROP,
};
PROP_OPT[] = {
    {"Noise", PROP_STR, {PROP_NO_VAL, "no"}, PROP_RNG_YESNO},
    {"Start", PROP_REAL, {1e9, PROP_NO_STR}, PROP_POS_RANGE},
    {"Stop", PROP_REAL, {10e9, PROP_NO_STR}, PROP_POS_RANGE},
    {"Points", PROP_INT, {10, PROP_NO_STR}, PROP_MIN_VAL(2)},
    {"Values", PROP_LIST, {10, PROP_NO_STR}, PROP_POS_RANGE},
    PROP_NO_PROP,
};
struct define_t acsolver::anadef = {"AC", 0, PROP_ACTION, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF};

} // namespace qucs
