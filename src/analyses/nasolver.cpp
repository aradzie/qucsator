/*
 * nasolver.cpp - nodal analysis solver class implementation
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

#include <cstring>

#include "analysis.h"
#include "circuit.h"
#include "complex.h"
#include "eqnsys.h"
#include "exception.h"
#include "exceptionstack.h"
#include "logging.h"
#include "nasolver.h"
#include "net.h"
#include "node.h"
#include "nodelist.h"
#include "nodeset.h"
#include "object.h"
#include "tmatrix.h"
#include "tvector.h"
#include "vector.h"

namespace qucs {

template <class nr_type_t> nasolver<nr_type_t>::nasolver() : analysis() {
  nlist = nullptr;
  A = C = nullptr;
  z = x = xprev = zprev = nullptr;
  reltol = abstol = vntol = 0;
  calculate_func = nullptr;
  convHelper = fixpoint = 0;
  eqnAlgo = ALGO_LU_DECOMPOSITION;
  updateMatrix = 1;
  gMin = srcFactor = 0;
  eqns = new eqnsys<nr_type_t>();
}

template <class nr_type_t> nasolver<nr_type_t>::nasolver(const std::string &n) : analysis(n) {
  nlist = nullptr;
  A = C = nullptr;
  z = x = xprev = zprev = nullptr;
  reltol = abstol = vntol = 0;
  calculate_func = nullptr;
  convHelper = fixpoint = 0;
  eqnAlgo = ALGO_LU_DECOMPOSITION;
  updateMatrix = 1;
  gMin = srcFactor = 0;
  eqns = new eqnsys<nr_type_t>();
}

template <class nr_type_t> nasolver<nr_type_t>::~nasolver() {
  delete nlist;
  delete C;
  delete A;
  delete z;
  delete x;
  delete xprev;
  delete zprev;
  delete eqns;
}

/* Runs the nodal analysis solver once, reports errors if any
 * and saves the results into each circuit. */
template <class nr_type_t> int nasolver<nr_type_t>::solve_once() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_once()\n", getName());

  int error = 0;

  // run the calculation function such as calcDC, calcTR, etc, for each circuit
  calculate();

  // generate A matrix and z vector
  createMatrix();

  // solve equation system
  runMNA();

  if (estack.top()) {
    estack.print();
    error = -1;
  }

  // save results into circuits
  if (!error) {
    saveSolution();
  }

  return error;
}

/* Run this function before the actual solver. */
template <class nr_type_t> void nasolver<nr_type_t>::solve_pre() {
  // create node list, enumerate nodes and voltage sources
#if DEBUG
  logprint(LOG_STATUS, "NOTIFY: %s: creating node list for %s analysis\n", getName(), desc.c_str());
#endif
  nlist = new nodelist(subnet);
  nlist->assignNodes();
  assignVoltageSources();
#if DEBUG
  nlist->print();
#endif

  // create matrix, solution vector and right hand side vector
  const int M = countVoltageSources();
  const int N = countNodes();
  delete A;
  A = new tmatrix<nr_type_t>(M + N);
  delete z;
  z = new tvector<nr_type_t>(N + M);
  delete x;
  x = new tvector<nr_type_t>(N + M);

#if DEBUG
  logprint(LOG_STATUS, "NOTIFY: %s: solving %s netlist\n", getName(), desc.c_str());
#endif
}

/* Run this function after the actual solver run and before evaluating
   the results. */
template <class nr_type_t> void nasolver<nr_type_t>::solve_post() {
  delete nlist;
  nlist = nullptr;
}

/* Goes through the nodeset list of the current netlist and applies the stored values
 * to the current solution vector.
 * Then the function saves the solution vector back into the actual component nodes. */
template <class nr_type_t> void nasolver<nr_type_t>::applyNodeset(bool nokeep) {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::applyNodeset()\n", getName());

  if (x == nullptr || nlist == nullptr)
    return;

  // set each solution to zero
  if (nokeep) {
    for (int i = 0; i < x->size(); i++) {
      x->set(i, 0);
    }
  }

  // then apply the nodeset itself
  for (nodeset *n = subnet->getNodeset(); n; n = n->getNext()) {
    struct nodelist_t *nl = nlist->getNode(n->getName());
    if (nl != nullptr) {
      x->set(nl->n, n->getValue());
    } else {
      logprint(LOG_ERROR,
               "WARNING: %s: no such node `%s' found, cannot "
               "initialize node\n",
               getName(), n->getName());
    }
  }
  if (xprev != nullptr) {
    *xprev = *x;
  }
  saveSolution();

  // propagate the solution to the non-linear circuits
  restartNR();
}

/* Uses the gMin-stepping algorithm in order to solve the given non-linear netlist
 * by continuous iterations. */
template <class nr_type_t> int nasolver<nr_type_t>::solve_nonlinear_continuation_gMin() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_nonlinear_continuation_gMin()\n", getName());

  const int MaxIter = getPropertyInteger("MaxIter") / 4 + 1;

  int error = 0;

  updateMatrix = 1;
  fixpoint = 0;

  // initialize the stepper
  double gPrev = gMin = 0.01;
  double gStep = gMin / 100;
  gMin -= gStep;

  do {
    // run solving loop until convergence is reached
    int iter = 0;
    bool convergence = false;
    do {
      error = solve_once();
      if (error) {
        break;
      }
      // convergence check
      convergence = (iter > 0) ? checkConvergence() : false;
      savePreviousIteration();
      iter++;
    } while (!convergence && iter < MaxIter);
    iterations += iter;

    // not yet converged, so decreased the gMin-step
    if (iter >= MaxIter || error) {
      gStep /= 2;
      // here the absolute minimum step checker
      if (gStep < std::numeric_limits<double>::epsilon()) {
        error = 1;
        const auto e = new qucs::exception(EXCEPTION_NO_CONVERGENCE);
        e->setText("no convergence in %s analysis after %d gMinStepping "
                   "iterations",
                   desc.c_str(), iterations);
        estack.push(e);
        break;
      }
      gMin = MAX(gPrev - gStep, 0);
    }
    // converged, increased the gMin-step
    else {
      gPrev = gMin;
      gMin = MAX(gMin - gStep, 0);
      gStep *= 2;
    }
  } while (gPrev > 0); // continue until no additional resistances is necessary

  return error;
}

/* The following function uses the source-stepping algorithm in order
   to solve the given non-linear netlist by continuous iterations. */
template <class nr_type_t> int nasolver<nr_type_t>::solve_nonlinear_continuation_Source() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_nonlinear_continuation_Source()\n", getName());

  const int MaxIter = getPropertyInteger("MaxIter") / 4 + 1;

  int error = 0;

  updateMatrix = 1;
  fixpoint = 0;

  // initialize the stepper
  double sPrev = srcFactor = 0;
  double sStep = 0.01;
  srcFactor += sStep;

  do {
    // run solving loop until convergence is reached
    int iter = 0;
    bool convergence;
    do {
      subnet->setSrcFactor(srcFactor);
      error = solve_once();
      if (error) {
        break;
      }
      // convergence check
      convergence = (iter > 0) ? checkConvergence() : false;
      savePreviousIteration();
      iter++;
    } while (!convergence && iter < MaxIter);
    iterations += iter;

    // not yet converged, so decreased the source-step
    if (iter >= MaxIter || error) {
      if (error) {
        sStep *= 0.1;
      } else {
        sStep *= 0.5;
      }
      restorePreviousIteration();
      saveSolution();
      // here the absolute minimum step checker
      if (sStep < std::numeric_limits<double>::epsilon()) {
        error = 1;
        const auto e = new qucs::exception(EXCEPTION_NO_CONVERGENCE);
        e->setText("no convergence in %s analysis after %d sourceStepping "
                   "iterations",
                   desc.c_str(), iterations);
        estack.push(e);
        break;
      }
      srcFactor = std::min(sPrev + sStep, 1.0);
    }
    // converged, increased the source-step
    else if (iter < MaxIter / 4) {
      sPrev = srcFactor;
      srcFactor = std::min(srcFactor + sStep, 1.0);
      sStep *= 1.5;
    } else {
      srcFactor = std::min(srcFactor + sStep, 1.0);
    }
  } while (sPrev < 1); // continue until no source factor is necessary

  subnet->setSrcFactor(1);
  return error;
}

/* The non-linear iterative nodal analysis netlist solver. */
template <class nr_type_t> int nasolver<nr_type_t>::solve_nonlinear() {
  const int MaxIter = getPropertyInteger("MaxIter");
  reltol = getPropertyDouble("reltol");
  abstol = getPropertyDouble("abstol");
  vntol = getPropertyDouble("vntol");

  int error = 0;

  updateMatrix = 1;

  if (convHelper == CONV_GMinStepping) {
    // use the alternative non-linear solver solve_nonlinear_continuation_gMin
    // instead of the basic solver provided by this function
    iterations = 0;
    error = solve_nonlinear_continuation_gMin();
    return error;
  }

  if (convHelper == CONV_SourceStepping) {
    // use the alternative non-linear solver solve_nonlinear_continuation_Source
    // instead of the basic solver provided by this function
    iterations = 0;
    error = solve_nonlinear_continuation_Source();
    return error;
  }

  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_nonlinear()\n", getName());

  log_indent();

  // run solving loop until convergence is reached
  int iter = 0;
  bool convergence;
  do {
    logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_nonlinear(), iter=%d\n", getName(), iter);

    error = solve_once();
    if (error) {
      break;
    }

    // convergence check
    convergence = (iter > 0) ? checkConvergence() : false;
    savePreviousIteration();
    iter++;
    // control fixpoint iterations
    if (fixpoint) {
      if (convergence && !updateMatrix) {
        updateMatrix = 1;
        convergence = false;
      } else {
        updateMatrix = 0;
      }
    }
  } while (!convergence && iter < MaxIter * (1 + convHelper ? 1 : 0));

  if (iter >= MaxIter || error) {
    const auto e = new qucs::exception(EXCEPTION_NO_CONVERGENCE);
    e->setText("no convergence in %s analysis after %d iterations", desc.c_str(), iter);
    estack.push(e);
    error++;
  }

  log_dedent();

  iterations = iter;
  return error;
}

/* The linear nodal analysis netlist solver. */
template <class nr_type_t> int nasolver<nr_type_t>::solve_linear() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::solve_linear(), iter=%d\n", getName());

  updateMatrix = 1;
  return solve_once();
}

/* Applying the MNA (Modified Nodal Analysis) to a circuit with
   passive elements and independent current and voltage sources
   results in a matrix equation of the form Ax = z.  This function
   generates the A and z matrix. */
template <class nr_type_t> void nasolver<nr_type_t>::createMatrix() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::createMatrix()\n", getName());

  /* Generate the A matrix.  The A matrix consists of four (4) minor
     matrices in the form     +-   -+
                          A = | G B |
                              | C D |
                              +-   -+.
     Each of these minor matrices is going to be generated here. */
  if (updateMatrix) {
    createGMatrix();
    createBMatrix();
    createCMatrix();
    createDMatrix();
  }

  /* Adjust G matrix if requested. */
  if (convHelper == CONV_GMinStepping) {
    const int N = countNodes();
    const int M = countVoltageSources();
    for (int n = 0; n < N + M; n++) {
      A->set(n, n, A->get(n, n) + gMin);
    }
  }

  /* Generate the z Matrix.  The z Matrix consists of two (2) minor
     matrices in the form     +- -+
                          z = | i |
                              | e |
                              +- -+.
     Each of these minor matrices is going to be generated here. */
  createZVector();
}

/* A helper to get the correct values from the circuit's matrices.
 * The additional (unused) argument is used to differentiate between the two possible types. */
#define MatVal(x) MatValX(x, (nr_type_t *)0)

template <class nr_type_t> nr_type_t nasolver<nr_type_t>::MatValX(nr_complex_t z, nr_complex_t *) {
  return z;
}

template <class nr_type_t> nr_type_t nasolver<nr_type_t>::MatValX(nr_complex_t z, double *) {
  return real(z);
}

/* The B matrix is an MxN matrix with only 0, 1 and -1 elements.  Each
   location in the matrix corresponds to a particular voltage source
   (first dimension) or a node (second dimension).  If the positive
   terminal of the ith voltage source is connected to node k, then the
   element (i,k) in the B matrix is a 1.  If the negative terminal of
   the ith voltage source is connected to node k, then the element
   (i,k) in the B matrix is a -1.  Otherwise, elements of the B matrix
   are zero. */
template <class nr_type_t> void nasolver<nr_type_t>::createBMatrix() {
  const int N = countNodes();
  const int M = countVoltageSources();

  // go through each voltage sources (first dimension)
  for (int c = 0; c < M; c++) {
    circuit *vs = findVoltageSource(c);
    // go through each node (second dimension)
    for (int r = 0; r < N; r++) {
      nr_type_t val = 0.0;
      struct nodelist_t *n = nlist->getNode(r);
      for (auto &current : *n) {
        // is voltage source connected to node ?
        if (current->getCircuit() == vs) {
          val += MatVal(vs->getB(current->getPort(), c));
        }
      }
      // put value into B matrix
      A->set(r, c + N, val);
    }
  }
}

/* The C matrix is an NxM matrix with only 0, 1 and -1 elements.  Each
   location in the matrix corresponds to a particular node (first
   dimension) or a voltage source (first dimension).  If the positive
   terminal of the ith voltage source is connected to node k, then the
   element (k,i) in the C matrix is a 1.  If the negative terminal of
   the ith voltage source is connected to node k, then the element
   (k,i) in the C matrix is a -1.  Otherwise, elements of the C matrix
   are zero. */
template <class nr_type_t> void nasolver<nr_type_t>::createCMatrix() {
  const int N = countNodes();
  const int M = countVoltageSources();
  // go through each voltage sources (second dimension)
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    // go through each node (first dimension)
    for (int c = 0; c < N; c++) {
      nr_type_t val = 0.0;
      struct nodelist_t *n = nlist->getNode(c);
      for (auto &current : *n) {
        // is voltage source connected to node ?
        if (current->getCircuit() == vs) {
          val += MatVal(vs->getC(r, current->getPort()));
        }
      }
      // put value into C matrix
      A->set(r + N, c, val);
    }
  }
}

/* The D matrix is an MxM matrix that is composed entirely of zeros.
   It can be non-zero if dependent sources are considered. */
template <class nr_type_t> void nasolver<nr_type_t>::createDMatrix() {
  const int M = countVoltageSources();
  const int N = countNodes();
  for (int r = 0; r < M; r++) {
    circuit *vsr = findVoltageSource(r);
    for (int c = 0; c < M; c++) {
      circuit *vsc = findVoltageSource(c);
      nr_type_t val = 0.0;
      if (vsr == vsc) {
        val = MatVal(vsr->getD(r, c));
      }
      A->set(r + N, c + N, val);
    }
  }
}

/* The G matrix is an NxN matrix formed in two steps.
   1. Each element in the diagonal matrix is equal to the sum of the
   conductance of each element connected to the corresponding node.
   2. The off diagonal elements are the negative conductance of the
   element connected to the pair of corresponding nodes.  Therefore a
   resistor between nodes 1 and 2 goes into the G matrix at location
   (1,2) and location (2,1).  If an element is grounded, it will only
   have contribute to one entry in the G matrix -- at the appropriate
   location on the diagonal. */
template <class nr_type_t> void nasolver<nr_type_t>::createGMatrix() {
  const int N = countNodes();
  // go through each column of the G matrix
  for (int c = 0; c < N; c++) {
    struct nodelist_t *nc = nlist->getNode(c);
    // go through each row of the G matrix
    for (int r = 0; r < N; r++) {
      struct nodelist_t *nr = nlist->getNode(r);
      nr_type_t g = 0.0;
      // sum up the conductance of each connected circuit
      for (auto &currentnc : *nc)
        for (auto &currentnr : *nr)
          if (currentnc->getCircuit() == currentnr->getCircuit()) {
            circuit *ct = currentnc->getCircuit();
            const int pc = currentnc->getPort();
            const int pr = currentnr->getPort();
            g += MatVal(ct->getY(pr, pc));
          }
      // put value into G matrix
      A->set(r, c, g);
    }
  }
}

/* The following function creates the (N+M)x(N+M) noise current
   correlation matrix used during the AC noise computations.  */
template <class nr_type_t> void nasolver<nr_type_t>::createNoiseMatrix() {
  const int N = countNodes();
  const int M = countVoltageSources();

  // create new Cy matrix if necessary
  delete C;
  C = new tmatrix<nr_type_t>(N + M);

  // go through each column of the Cy matrix
  for (int c = 0; c < N; c++) {
    struct nodelist_t *nc = nlist->getNode(c);
    // go through each row of the Cy matrix
    for (int r = 0; r < N; r++) {
      struct nodelist_t *nr = nlist->getNode(r);
      nr_type_t val = 0.0;
      // sum up the noise-correlation of each connected circuit
      for (auto &currentnc : *nc)
        /* a = 0; a < nc->size(); a++ */
        for (auto &currentnr : *nr)
          /* b = 0; b < nr->size(); b++) */
          if (currentnc->getCircuit() == currentnr->getCircuit()) {
            circuit *ct = currentnc->getCircuit();
            int pc = currentnc->getPort();
            int pr = currentnr->getPort();
            val += MatVal(ct->getN(pr, pc));
          }
      // put value into Cy matrix
      C->set(r, c, val);
    }
  }

  // go through each additional voltage source and put coefficients into
  // the noise current correlation matrix
  circuit *vsr, *vsc;
  for (int r = 0; r < M; r++) {
    vsr = findVoltageSource(r);
    for (int c = 0; c < M; c++) {
      vsc = findVoltageSource(c);
      nr_type_t val = 0.0;
      if (vsr == vsc) {
        int ri = vsr->getSize() + r - vsr->getVoltageSource();
        int ci = vsc->getSize() + c - vsc->getVoltageSource();
        val = MatVal(vsr->getN(ri, ci));
      }
      C->set(r + N, c + N, val);
    }
  }

  // go through each additional voltage source
  for (int r = 0; r < M; r++) {
    vsr = findVoltageSource(r);
    // go through each node
    for (int c = 0; c < N; c++) {
      nr_type_t val = 0.0;
      struct nodelist_t *n = nlist->getNode(c);
      for (auto &currentn : *n)
      /*i = 0; i < n->size(); i++ )*/
      {
        // is voltage source connected to node ?
        if (currentn->getCircuit() == vsr) {
          int ri = vsr->getSize() + r - vsr->getVoltageSource();
          int ci = currentn->getPort();
          val += MatVal(vsr->getN(ri, ci));
        }
      }
      // put value into Cy matrix
      C->set(r + N, c, val);
    }
  }

  // go through each voltage source
  for (int c = 0; c < M; c++) {
    vsc = findVoltageSource(c);
    // go through each node
    for (int r = 0; r < N; r++) {
      nr_type_t val = 0.0;
      struct nodelist_t *n = nlist->getNode(r);
      for (auto &currentn : *n) /*i = 0; i < n->size(); i++)*/
      {
        // is voltage source connected to node ?
        if (currentn->getCircuit() == vsc) {
          int ci = vsc->getSize() + c - vsc->getVoltageSource();
          int ri = currentn->getPort();
          val += MatVal(vsc->getN(ri, ci));
        }
      }
      // put value into Cy matrix
      C->set(r, c + N, val);
    }
  }
}

/* The i matrix is an 1xN matrix with each element of the matrix
   corresponding to a particular node.  The value of each element of i
   is determined by the sum of current sources into the corresponding
   node.  If there are no current sources connected to the node, the
   value is zero. */
template <class nr_type_t> void nasolver<nr_type_t>::createIVector() {
  const int N = countNodes();
  // go through each node
  for (int r = 0; r < N; r++) {
    nr_type_t val = 0.0;
    struct nodelist_t *n = nlist->getNode(r);
    // go through each circuit connected to the node
    for (auto &currentn : *n) /* int i = 0; i < n->size(); i++)*/
    {
      circuit *is = currentn->getCircuit();
      // is this a current source ?
      if (is->isISource() || is->isNonLinear()) {
        val += MatVal(is->getI(currentn->getPort()));
      }
    }
    // put value into i vector
    z->set(r, val);
  }
}

/* The e matrix is a 1xM matrix with each element of the matrix equal
   in value to the corresponding independent voltage source. */
template <class nr_type_t> void nasolver<nr_type_t>::createEVector() {
  const int N = countNodes();
  const int M = countVoltageSources();
  // go through each voltage source
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    nr_type_t val = MatVal(vs->getE(r));
    // put value into e vector
    z->set(r + N, val);
  }
}

// The function loads the right hand side vector.
template <class nr_type_t> void nasolver<nr_type_t>::createZVector() {
  createIVector();
  createEVector();
}

// Returns the number of nodes in the nodelist, excluding the ground node.
template <class nr_type_t> int nasolver<nr_type_t>::countNodes() { return nlist->length() - 1; }

// Returns the node number of the give node name.
template <class nr_type_t> int nasolver<nr_type_t>::getNodeNr(const std::string &str) {
  return nlist->getNodeNr(str);
}

/* Returns the assigned node number for the port of the given circuits,
 * or -1 if there is no such node. */
template <class nr_type_t> int nasolver<nr_type_t>::findAssignedNode(circuit *c, int port) {
  const int N = countNodes();
  for (int r = 0; r < N; r++) {
    struct nodelist_t *n = nlist->getNode(r);
    for (auto &currentn : *n) /*int i = 0; i < n->size(); i++)*/
      if (c == currentn->getCircuit())
        if (port == currentn->getPort())
          return r;
  }
  return -1;
}

// Returns the number of voltage sources in the nodelist.
template <class nr_type_t> int nasolver<nr_type_t>::countVoltageSources() {
  return subnet->getVoltageSources();
}

/* Returns the voltage source circuit object corresponding to the given number.
 * If there is no such voltage source it returns nullptr. */
template <class nr_type_t> circuit *nasolver<nr_type_t>::findVoltageSource(int n) {
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (n >= c->getVoltageSource() && n <= c->getVoltageSource() + c->getVoltageSources() - 1)
      return c;
  }
  return nullptr;
}

/* Applies unique voltage source identifiers to each voltage source
 * (explicit and built in internal ones) in the list of registered circuits. */
template <class nr_type_t> void nasolver<nr_type_t>::assignVoltageSources() {
  circuit *root = subnet->getRoot();
  int nSources = 0;
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->getVoltageSources() > 0) {
      c->setVoltageSource(nSources);
      nSources += c->getVoltageSources();
    }
  }
  subnet->setVoltageSources(nSources);
}

/* The matrix equation Ax = z is solved by x = A^-1*z.  The function
   applies the operation to the previously generated matrices. */
template <class nr_type_t> void nasolver<nr_type_t>::runMNA() {
  // just solve the equation system here
  eqns->setAlgo(eqnAlgo);
  eqns->passEquationSys(updateMatrix ? A : nullptr, x, z);
  eqns->solve();

  // if damped Newton-Raphson is requested
  if (xprev != nullptr && estack.top() == nullptr) {
    if (convHelper == CONV_Attenuation) {
      applyAttenuation();
    } else if (convHelper == CONV_LineSearch) {
      lineSearch();
    } else if (convHelper == CONV_SteepestDescent) {
      steepestDescent();
    }
  }
}

/* This function applies a damped Newton-Raphson (limiting scheme) to
   the current solution vector in the form x1 = x0 + a * (x1 - x0).  This
   convergence helper is heuristic and does not ensure global convergence. */
template <class nr_type_t> void nasolver<nr_type_t>::applyAttenuation() {
  double alpha = 1.0;

  // create solution difference vector and find maximum deviation
  tvector<nr_type_t> dx = *x - *xprev;
  double nMax = maxnorm(dx);

  // compute appropriate damping factor
  if (nMax > 0.0) {
    double g = 1.0;
    alpha = std::min(0.9, g / nMax);
    if (alpha < 0.1)
      alpha = 0.1;
  }

  // apply damped solution vector
  *x = *xprev + alpha * dx;
}

/* This is damped Newton-Raphson using nested iterations in order to
   find a better damping factor.  It identifies a damping factor in
   the interval [0,1] which minimizes the right hand side vector.  The
   algorithm actually ensures global convergence but pushes the
   solution to local minimums, i.e. where the Jacobian matrix A may be
   singular. */
template <class nr_type_t> void nasolver<nr_type_t>::lineSearch() {
  double alpha = 0.5, n, nMin, aprev = 1.0, astep = 0.5, adiff;
  int dir = -1;

  // compute solution deviation vector
  tvector<nr_type_t> dx = *x - *xprev;
  nMin = std::numeric_limits<double>::max();

  do {
    // apply current damping factor and see what happens
    *x = *xprev + alpha * dx;

    // recalculate Jacobian and right hand side
    saveSolution();
    calculate();
    createZVector();

    // calculate norm of right hand side vector
    n = norm(*z);

    // TODO: this is not perfect, but usable
    astep /= 2;
    adiff = fabs(alpha - aprev);
    if (adiff > 0.005) {
      aprev = alpha;
      if (n < nMin) {
        nMin = n;
        if (alpha == 1)
          dir = -dir;
        alpha += astep * dir;
      } else {
        dir = -dir;
        alpha += 1.5 * astep * dir;
      }
    }
  } while (adiff > 0.005);

  // apply final damping factor
  assert(alpha > 0 && alpha <= 1);
  *x = *xprev + alpha * dx;
}

/* The function looks for the optimal gradient for the right hand side
   vector using the so-called 'steepest descent' method.  Though
   better than the one-dimensional linesearch (it doesn't push
   iterations into local minimums) it converges painfully slow. */
template <class nr_type_t> void nasolver<nr_type_t>::steepestDescent() {
  double alpha = 1.0, sl, n;

  // compute solution deviation vector
  tvector<nr_type_t> dx = *x - *xprev;
  tvector<nr_type_t> dz = *z - *zprev;
  n = norm(*zprev);

  do {
    // apply current damping factor and see what happens
    *x = *xprev + alpha * dx;

    // recalculate Jacobian and right hand side
    saveSolution();
    calculate();
    createZVector();

    // check gradient criteria, ThinkME: Is this correct?
    dz = *z - *zprev;
    sl = real(sum(dz * -dz));
    if (norm(*z) < n + alpha * sl)
      break;
    alpha *= 0.7;
  } while (alpha > 0.001);

  // apply final damping factor
  *x = *xprev + alpha * dx;
}

/* Checks whether the iterative algorithm for linearizing the non-linear components
 * in the network shows convergence.
 * It returns non-zero if it converges and zero otherwise. */
template <class nr_type_t> bool nasolver<nr_type_t>::checkConvergence() {
  const int N = countNodes();
  const int M = countVoltageSources();

  // check the nodal voltage changes against the allowed absolute and relative tolerance values
  for (int r = 0; r < N; r++) {
    const double v_abs = abs(x->get(r) - xprev->get(r));
    const double v_rel = abs(x->get(r));
    if (v_abs >= vntol + reltol * v_rel) {
      return false;
    }
    if (!convHelper) {
      const double i_abs = abs(z->get(r) - zprev->get(r));
      const double i_rel = abs(z->get(r));
      if (i_abs >= abstol + reltol * i_rel) {
        return false;
      }
    }
  }

  for (int r = 0; r < M; r++) {
    const double i_abs = abs(x->get(r + N) - xprev->get(r + N));
    const double i_rel = abs(x->get(r + N));
    if (i_abs >= abstol + reltol * i_rel) {
      return false;
    }
    if (!convHelper) {
      const double v_abs = abs(z->get(r + N) - zprev->get(r + N));
      const double v_rel = abs(z->get(r + N));
      if (v_abs >= vntol + reltol * v_rel) {
        return false;
      }
    }
  }

  return true;
}

/* Saves the solution and right hand vector of the previous iteration. */
template <class nr_type_t> void nasolver<nr_type_t>::savePreviousIteration() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::savePreviousIteration()\n", getName());

  if (xprev != nullptr) {
    *xprev = *x;
  } else {
    xprev = new tvector<nr_type_t>(*x);
  }
  if (zprev != nullptr) {
    *zprev = *z;
  } else {
    zprev = new tvector<nr_type_t>(*z);
  }
}

/* Restores the solution and right hand vector of the previous (successful) iteration. */
template <class nr_type_t> void nasolver<nr_type_t>::restorePreviousIteration() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::restorePreviousIteration()\n", getName());

  if (xprev != nullptr) {
    *x = *xprev;
  }
  if (zprev != nullptr) {
    *z = *zprev;
  }
}

/* Restarts the NR iteration for each non-linear circuit. */
template <class nr_type_t> void nasolver<nr_type_t>::restartNR() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::restartDC()\n", getName());

  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->isNonLinear()) {
      c->restartDC();
    }
  }
}

/* Goes through solution (the x vector) and saves the node voltages of the last iteration
 * into each non-linear circuit. */
template <class nr_type_t> void nasolver<nr_type_t>::saveNodeVoltages() {
  const int N = countNodes();
  struct nodelist_t *n;
  // save all nodes except reference node
  for (int r = 0; r < N; r++) {
    n = nlist->getNode(r);
    /* for (int i = 0; i < n->size(); i++)*/
    for (auto &currentn : *n) {
      currentn->getCircuit()->setV(currentn->getPort(), x->get(r));
    }
  }
  // save reference node
  n = nlist->getNode(-1);
  for (auto &currentn : *n)
    currentn->getCircuit()->setV(currentn->getPort(), 0.0);
}

/* Goes through solution (the x vector) and saves the branch currents through the voltage sources
 * of the last iteration into each circuit. */
template <class nr_type_t> void nasolver<nr_type_t>::saveBranchCurrents() {
  const int N = countNodes();
  const int M = countVoltageSources();
  circuit *vs;
  // save all branch currents of voltage sources
  for (int r = 0; r < M; r++) {
    vs = findVoltageSource(r);
    vs->setJ(r, x->get(r + N));
  }
}

// The function saves the solution vector into each circuit.
template <class nr_type_t> void nasolver<nr_type_t>::saveSolution() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::saveSolution()\n", getName());

  saveNodeVoltages();
  saveBranchCurrents();
}

// This function stores the solution (node voltages and branch currents).
template <class nr_type_t> void nasolver<nr_type_t>::storeSolution() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::storeSolution()\n", getName());

  // cleanup solution previously
  solution.clear();
  const int N = countNodes();
  const int M = countVoltageSources();
  // store all nodes except reference node
  for (int r = 0; r < N; r++) {
    struct nodelist_t *n = nlist->getNode(r);
    nr_type_t gr = x->get(r);
    naentry<nr_type_t> entry(gr, 0);
    solution.insert({{n->name, entry}});
  }
  // store all branch currents of voltage sources
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    int vn = r - vs->getVoltageSource() + 1;
    nr_type_t xg = x->get(r + N);
    naentry<nr_type_t> entry(xg, vn);
    solution.insert({{vs->getName(), entry}});
  }
}

// This function recalls the solution (node voltages and branch currents).
template <class nr_type_t> void nasolver<nr_type_t>::recallSolution() {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::recallSolution()\n", getName());

  const int N = countNodes();
  const int M = countVoltageSources();
  // store all nodes except reference node
  for (int r = 0; r < N; r++) {
    struct nodelist_t *n = nlist->getNode(r);
    auto na = solution.find(n->name);
    if (na != solution.end())
      if ((*na).second.current == 0)
        x->set(r, (*na).second.value);
  }
  // store all branch currents of voltage sources
  for (int r = 0; r < M; r++) {
    circuit *vs = findVoltageSource(r);
    int vn = r - vs->getVoltageSource() + 1;
    auto na = solution.find(vs->getName());
    if (na != solution.end())
      if ((*na).second.current == vn)
        x->set(r + N, (*na).second.value);
  }
}

/* Saves the results of a single solve() functionality into the output dataset. */
template <class nr_type_t>
void nasolver<nr_type_t>::saveResults(const std::string &volts, const std::string &amps,
                                      int saveOPs, qucs::vector *f) {
  logprint(LOG_STATUS, "NOTIFY: %s: nasolver::saveResults()\n", getName());

  const int N = countNodes();
  const int M = countVoltageSources();

  // add node voltage variables
  if (!volts.empty()) {
    for (int r = 0; r < N; r++) {
      std::string n = createV(r, volts, saveOPs);
      if (!n.empty()) {
        saveVariable(n, x->get(r), f);
      }
    }
  }

  // add branch current variables
  if (!amps.empty()) {
    for (int r = 0; r < M; r++) {
      std::string n = createI(r, amps, saveOPs);
      if (!n.empty()) {
        saveVariable(n, x->get(r + N), f);
      }
    }
  }

  // add voltage probe data
  if (!volts.empty()) {
    circuit *root = subnet->getRoot();
    for (circuit *c = root; c != nullptr; c = c->getNext()) {
      if (!c->isProbe())
        continue;
      if (!c->getSubcircuit().empty() && !(saveOPs & SAVE_ALL))
        continue;
      if (volts != "vn")
        c->saveOperatingPoints();
      std::string n = createOP(c->getName(), volts);
      saveVariable(n, nr_complex_t(c->getOperatingPoint("Vr"), c->getOperatingPoint("Vi")), f);

      // add watt probe data
      c->calcOperatingPoints();
      for (auto ops : c->getOperatingPoints()) {
        // It will only get values if none of the strings are 0
        // Once again most of this is adapted from Vprobe and Iprobe
        qucs::pair &p = ops.second;
        if (strcmp(p.getName(), "Vi") == 0)
          continue;
        if (strcmp(p.getName(), "VAi") == 0)
          continue;
        if (strcmp(p.getName(), "Vr") == 0)
          continue;
        if (strcmp(p.getName(), "VAr") == 0) {
          std::string n = createOP(c->getName(), "S");
          saveVariable(n, nr_complex_t(c->getOperatingPoint("VAr"), c->getOperatingPoint("VAi")),
                       f);
          continue;
        }

        std::string n = createOP(c->getName(), p.getName());
        saveVariable(n, p.getValue(), f);
      }
    }
  }

  // save operating points of non-linear circuits if requested
  if (saveOPs & SAVE_OPS) {
    circuit *root = subnet->getRoot();
    for (circuit *c = root; c != nullptr; c = c->getNext()) {
      if (!c->isNonLinear())
        continue;
      if (!c->getSubcircuit().empty() && !(saveOPs & SAVE_ALL))
        continue;
      c->calcOperatingPoints();
      for (auto ops : c->getOperatingPoints()) {
        qucs::pair &p = ops.second;
        std::string n = createOP(c->getName(), p.getName());
        saveVariable(n, p.getValue(), f);
      }
    }
  }
}

/* Creates an appropriate variable name for operating points.
 * The caller is responsible to free() the returned string. */
template <class nr_type_t>
std::string nasolver<nr_type_t>::createOP(const std::string &c, const std::string &n) {
  return c + "." + n;
}

/* Creates an appropriate variable name for voltages.
 * The caller is responsible to free() the returned string. */
template <class nr_type_t>
std::string nasolver<nr_type_t>::createV(int n, const std::string &volts, int saveOPs) {
  if (nlist->isInternal(n))
    return std::string();
  std::string node = nlist->get(n);
  if (node.find('.') != std::string::npos && !(saveOPs & SAVE_ALL))
    return std::string();
  std::string ret = node + "." + volts;
  return ret;
}

/* Creates an appropriate variable name for currents.
 * The caller is responsible to free() the returned string. */
template <class nr_type_t>
std::string nasolver<nr_type_t>::createI(int n, const std::string &amps, int saveOPs) {
  circuit *vs = findVoltageSource(n);

  // don't output internal (helper) voltage sources
  if (vs->isInternalVoltageSource())
    return std::string();

  /* save only current through real voltage sources and explicit
     current probes */
  if (!vs->isVSource() && !(saveOPs & SAVE_OPS))
    return std::string();

  // don't output subcircuit components if not requested
  if (!vs->getSubcircuit().empty() && !(saveOPs & SAVE_ALL))
    return std::string();

  // create appropriate current name for single/multiple voltage sources
  std::string name = vs->getName();
  if (vs->getVoltageSources() > 1)
    return name + "." + amps + std::to_string(n - vs->getVoltageSource() + 1);
  else
    return name + "." + amps;
}

/* Alternative to countNodes () */
template <class nr_type_t> int nasolver<nr_type_t>::getN() { return countNodes(); }

/* Alternative to countVoltageSources () */
template <class nr_type_t> int nasolver<nr_type_t>::getM() { return countVoltageSources(); }

/* Returns an appropriate text representation for the currently used
 * convergence helper algorithm. */
template <class nr_type_t> const char *nasolver<nr_type_t>::getHelperDescription() {
  if (convHelper == CONV_Attenuation) {
    return "RHS attenuation";
  }
  if (convHelper == CONV_LineSearch) {
    return "line search";
  }
  if (convHelper == CONV_SteepestDescent) {
    return "steepest descent";
  }
  if (convHelper == CONV_GMinStepping) {
    return "gMin stepping";
  }
  if (convHelper == CONV_SourceStepping) {
    return "source stepping";
  }
  return "none";
}

} // namespace qucs
