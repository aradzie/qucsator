/*
 * net.cpp - net class implementation
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2007 Stefan Jahn <stefan@lkcc.org>
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

#include <cassert>
#include <cstring>

#include "analysis.h"
#include "circuit.h"
#include "component_id.h"
#include "dataset.h"
#include "environment.h"
#include "logging.h"
#include "net.h"
#include "node.h"
#include "nodelist.h"
#include "nodeset.h"
#include "object.h"
#include "ptrlist.h"

namespace qucs {

net::net() : object() {
  root = drop = nullptr;
  nPorts = nCircuits = nSources = 0;
  insertedNodes = inserted = reduced = 0;
  actions = new ptrlist<analysis>();
  orgacts = new ptrlist<analysis>();
  env = nullptr;
  nset = nullptr;
  srcFactor = 1;
}

net::net(const std::string &n) : object(n) {
  root = drop = nullptr;
  nPorts = nCircuits = nSources = 0;
  insertedNodes = inserted = reduced = 0;
  actions = new ptrlist<analysis>();
  orgacts = new ptrlist<analysis>();
  env = nullptr;
  nset = nullptr;
  srcFactor = 1;
}

net::~net() {
  circuit *n;
  // delete each and every circuit
  for (circuit *c = root; c != nullptr; c = n) {
    n = c->getNext();
    delete c;
  }
  // delete original actions
  for (auto *element : *orgacts) {
    delete element;
    element = nullptr;
  }

  delete orgacts;
  // delete nodeset
  delNodeset();
  delete actions;
}

/* Prepends the given circuit to the list of registered circuits. */
void net::insertCircuit(circuit *c) {
  assert(!containsCircuit(c));

  // chain circuit appropriately
  if (root)
    root->setPrev(c);
  c->setNext(root);
  c->setPrev(nullptr);
  root = c;
  nCircuits++;
  c->setEnabled(1);
  c->setNet(this);

  /* handle AC power sources as s-parameter ports if it is not part of
     a subcircuit */
  if (c->getType() == CIR_PAC && c->getSubcircuit().empty()) {
    nPorts++;
    if (!c->getPort())
      c->setPort(c->getPropertyInteger("Num"));
  }
  // handle DC voltage sources
  if (c->getVoltageSources() > 0) {
    if (c->getVoltageSource() < 0)
      c->setVoltageSource(nSources);
    nSources += c->getVoltageSources();
  }
}

/* Removes the given circuit from the list of registered circuits. */
void net::removeCircuit(circuit *c, int dropping) {
  assert(containsCircuit(c));

  // adjust the circuit chain appropriately
  if (c == root) {
    root = c->getNext();
    if (root)
      root->setPrev(nullptr);
  } else {
    if (c->getNext())
      c->getNext()->setPrev(c->getPrev());
    c->getPrev()->setNext(c->getNext());
  }
  nCircuits--;
  c->setEnabled(0);
  c->setNet(nullptr);
  if (c->getPort())
    nPorts--;
  if (c->getVoltageSource() >= 0)
    nSources -= c->getVoltageSources();

  // shift the circuit object to the drop list
  if (c->isOriginal()) {
    if (dropping) {
      if (drop)
        drop->setPrev(c);
      c->setNext(drop);
      c->setPrev(nullptr);
      drop = c;
    }
  } else {
    // really destroy the circuit object
    delete c;
  }
}

/* Returns non-zero if the given circuit is already part
   of the netlist. It returns zero if not. */
int net::containsCircuit(circuit *cand) {
  for (circuit *c = root; c != nullptr; c = c->getNext())
    if (c == cand)
      return 1;
  return 0;
}

/* This function prepends the given analysis to the list of registered
   analyses. */
void net::insertAnalysis(analysis *a) {
  orgacts->push_front(a);
  actions->push_front(a);
}

/* The function removes the given analysis from the list of registered
   analyses. */
void net::removeAnalysis(analysis *a) { actions->remove(a); }

/* Returns the analysis associated with the netlist
   object specified by the given instance name and returns nullptr if
   there is no such analysis. */
analysis *net::findAnalysis(const std::string &n) const {
  for (auto *a : *actions) {
    if (a->getName() == n)
      return a;
  }
  return nullptr;
}

/* Returns the analysis associated with the netlist
   object specified by the given type of analysis and returns nullptr if
   there is no such analysis. */
analysis *net::findAnalysis(int type) {
  for (auto *a : *actions) {
    if (a->getType() == type)
      return a;
  }
  return nullptr;
}

/* Looks recursively for a type of analysis. */
int net::containsAnalysis(analysis *child, int type) {
  ptrlist<analysis> *alist = child->getAnalysis();
  if (alist != nullptr) {
    for (auto *a : *alist) {
      if (a->getType() == type)
        return 1;
      else if (a->getType() == ANALYSIS_SWEEP)
        return containsAnalysis(a, type);
    }
  }
  return 0;
}

/* Runs all registered analyses applied to the current
   netlist, except for external analysis types. */
dataset *net::runAnalysis(int &err) {
  dataset *out = new dataset();

  // apply some data to all analyses
  for (auto *a : *actions) {
    a->setNet(this);
    a->setData(out);
  }

  // re-order analyses
  orderAnalysis();

  // initialize analyses
  for (auto *a : *actions) {
    err |= a->initialize();
  }

  // solve the analyses
  for (auto *a : *actions) {
    a->getEnv()->runSolver();
    err |= a->solve();
  }

  // cleanup analyses
  for (auto *a : *actions) {
    err |= a->cleanup();
  }

  return out;
}

/* Returns the analysis with the second lowest order.  If
   there is no recursive sweep it returns nullptr. */
analysis *net::findSecondOrder() {
  analysis *parent = nullptr;
  for (auto *a : *actions) {
    // parameter sweeps are potential parent sweeps
    if (a->getType() == ANALYSIS_SWEEP) {
      // find the appropriate sub analysis
      analysis *child = getChildAnalysis(a);
      if (child != nullptr) {
        // check if child is not another variable sweep
        if (child->getType() != ANALYSIS_SWEEP) {
          parent = a;
          break;
        }
        // check if the child's child is still in the analysis list
        else if (getChildAnalysis(child) == nullptr) {
          parent = a;
          break;
        }
      }
    }
  }
  return parent;
}

/* Reorders (prioritizes) the registered analysis to the netlist object.
 * It chains the analyses to be executed in a certain order. */
void net::orderAnalysis() {
  analysis *parent, *child;
  analysis *dc = findAnalysis(ANALYSIS_DC);
  int dcApplied = 0;
  do {
    // get second order sweep
    if ((parent = findSecondOrder()) != nullptr) {
      child = getChildAnalysis(parent);
      removeAnalysis(child);
      // apply sub-analysis to each parent analysis if any
      if (actions != nullptr) {
        for (auto *a : *actions) {
          const char *cn = getChild(a);
          if (cn != nullptr && !strcmp(cn, child->getName())) {
            a->addAnalysis(child);
            // apply DC analysis if necessary
            if (child->getType() != ANALYSIS_DC && child->getType() != ANALYSIS_SWEEP &&
                dc != nullptr) {
              if (!dcApplied)
                removeAnalysis(dc);
              a->addAnalysis(dc);
              dcApplied++;
            }
          }
        }
      }
      // sort the sub-analysis of each parent
      for (auto *a : *actions) {
        sortChildAnalyses(a);
      }
    }
  } while (parent != nullptr);

  // sort the parent analyses
  parent = new analysis();
  parent->setAnalysis(actions);
  sortChildAnalyses(parent);
  actions = new ptrlist<analysis>(*(parent->getAnalysis()));
  delete parent;
}

// Sorts the analyses of the given parent analysis.
void net::sortChildAnalyses(analysis *parent) {
  ptrlist<analysis> *alist = parent->getAnalysis();
  if (alist != nullptr) {

    for (auto it = alist->begin(); it != alist->end(); /* empty */) {
      // Copy the value of the element (a pointer), and advance the
      // iterator prior to manipulating the list.
      analysis *a = *it;
      ++it;

      if (a->getType() == ANALYSIS_DC || containsAnalysis(a, ANALYSIS_DC)) {
        parent->delAnalysis(a);
        parent->addAnalysis(a);
      }
    }
  }
}

// Returns the instance name of the given parents child analysis.
const char *net::getChild(analysis *parent) const {
  const char *child = nullptr;
  if (parent != nullptr && parent->getType() == ANALYSIS_SWEEP)
    child = parent->getPropertyString("Sim");
  return child;
}

// Returns the child analysis of the given parent if possible.
analysis *net::getChildAnalysis(analysis *parent) { return findAnalysis(getChild(parent)); }

// Returns the last order sweep being not an parameter sweep.
analysis *net::findLastOrder(analysis *a) {
  ptrlist<analysis> *alist = a->getAnalysis();
  analysis *child = alist ? alist->front() : nullptr;
  if (child != nullptr && child->getType() == ANALYSIS_SWEEP) {
    return findLastOrder(child);
  }
  return child ? child : a;
}

// Returns the last order sweep being not an parameter sweep.
ptrlist<analysis> *net::findLastOrderChildren(analysis *a) {
  ptrlist<analysis> *alist = a->getAnalysis();
  analysis *child = alist ? alist->front() : nullptr;
  if (child != nullptr && child->getType() == ANALYSIS_SWEEP) {
    return findLastOrderChildren(child);
  }
  return alist;
}

/* Re-shifts all circuits in the drop list to the actual list of circuit objects. */
void net::getDroppedCircuits(nodelist *nodes) {
  circuit *n;
  for (circuit *c = drop; c != nullptr; c = n) {
    n = c->getNext();
    if (nodes)
      nodes->insert(c);
    insertCircuit(c);
  }
  drop = nullptr;
}

/* Deletes all unnecessary circuits in the list of registered circuit objects. */
void net::deleteUnusedCircuits(nodelist *nodes) {
  circuit *n;
  for (circuit *c = root; c != nullptr; c = n) {
    n = c->getNext();
    if (!c->isOriginal()) {
      if (nodes)
        nodes->remove(c);
      removeCircuit(c);
    }
  }
}

/* Returns the first node in the list of real circuit objects
   connected to the given node.  If there is no such node (unconnected
   node) the function returns nullptr. */
node *net::findConnectedCircuitNode(node *n) {

  const char *_name = n->getName();
  node *_node;

  // through the list of circuit objects
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    // skip signal circuits
    if (c->getPort())
      continue;
    // through the list of nodes in a circuit
    for (int i = 0; i < c->getSize(); i++) {
      _node = c->getNode(i);
      if (!strcmp(_node->getName(), _name)) {
        if (_node != n) {
          return _node;
        }
      }
    }
  }
  return nullptr;
}

/* Returns the first node in the list of circuit objects (including
   signals) connected to the given node.  If there is no such node
   (unconnected node) the function returns nullptr. */
node *net::findConnectedNode(node *n) {

  const char *_name = n->getName();
  node *_node;

  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    for (int i = 0; i < c->getSize(); i++) {
      _node = c->getNode(i);
      if (!strcmp(_node->getName(), _name)) {
        if (_node != n) {
          return _node;
        }
      }
    }
  }
  return nullptr;
}

// Renames the given circuit and mark it as being a reduced one.
void net::reducedCircuit(circuit *c) {
  char n[32];
  sprintf(n, "reduced%d", reduced++);
  c->setName(n);
}

/* Renames the given circuit and mark it as being a inserted one and
   remember when it was inserted. */
void net::insertedCircuit(circuit *c) {
  char n[32];
  sprintf(n, "inserted%d", inserted);
  c->setName(n);
  c->setInserted(inserted);
  inserted++;
}

// Renames the given node and mark it as being a inserted one.
void net::insertedNode(node *c) {
  char n[32];
  sprintf(n, "inode%d", insertedNodes++);
  c->setName(n);
}

/* Checks whether the circuit chain of the netlist is properly working.
 * It returns the number of errors or zero if there are no errors. */
int net::checkCircuitChain() {
  int error = 0;
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->getPrev())
      if (c->getPrev()->getNext() != c) {
        error++;
        logprint(LOG_ERROR, "ERROR: prev->next != circuit '%s'\n", c->getName());
      }
    if (c->getNext())
      if (c->getNext()->getPrev() != c) {
        error++;
        logprint(LOG_ERROR, "ERROR: next->prev != circuit '%s'\n", c->getName());
      }
  }
  return error;
}

/* Counts the number of signals (ports) within the list of registered circuits. */
int net::countPorts() {
  int count = 0;
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->getPort()) {
      count++;
    }
  }
  return count;
}

/* Counts the number of circuits within the list of registered circuits. */
int net::countNodes() {
  int count = 0;
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (!c->getPort()) {
      count += c->getSize();
    }
  }
  return count;
}

/* Returns the number of non-linear circuits within the list of registered circuits. */
int net::isNonLinear() {
  int count = 0;
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->isNonLinear()) {
      count++;
    }
  }
  return count;
}

/* Adds the given nodeset object to the netlist's nodeset list. */
void net::addNodeset(nodeset *n) {
  n->setNext(nset);
  nset = n;
}

/* Deletes all the nodeset list of the netlist object.  Called from the destructor. */
void net::delNodeset() {
  nodeset *next;
  for (nodeset *n = nset; n != nullptr; n = next) {
    next = n->getNext();
    delete n;
  }
  nset = nullptr;
}

void net::setActionNetAll(net *subnet) {
  for (auto *a : *(this->actions)) {
    a->setNet(subnet);
  }
}

#if DEBUG
void net::list() {
  logprint(LOG_STATUS,
           "DEBUG: netlist `%s' (%d circuits, "
           "%d ports, %d nodes)\n",
           getName(), countPorts(), countPorts(), countNodes());
  // go through circuit list
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    // list each circuit
    logprint(LOG_STATUS, "       %s[", c->getName());
    for (int i = 0; i < c->getSize(); i++) {
      logprint(LOG_STATUS, "%s-%d", c->getNode(i)->getName(), c->getNode(i)->getNode());
      if (i < c->getSize() - 1)
        logprint(LOG_STATUS, ",");
    }
    logprint(LOG_STATUS, "] { %s }\n", c->propertyList());
  }
}
#endif /* DEBUG */

} // namespace qucs
