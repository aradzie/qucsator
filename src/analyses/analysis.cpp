/*
 * analysis.cpp - analysis class implementation
 *
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
#include "complex.h"
#include "dataset.h"
#include "object.h"
#include "ptrlist.h"
#include "strlist.h"
#include "sweep.h"
#include "vector.h"

namespace qucs {

analysis::analysis() : object() {
  data = nullptr;
  subnet = nullptr;
  env = nullptr;
  actions = nullptr;
  type = ANALYSIS_UNKNOWN;
  runs = 0;
  progress = true;
}

analysis::analysis(const std::string &n) : object(n) {
  data = nullptr;
  subnet = nullptr;
  env = nullptr;
  actions = nullptr;
  type = ANALYSIS_UNKNOWN;
  runs = 0;
  progress = true;
}

analysis::~analysis() { delete actions; }

analysis::analysis(analysis &a) : object(a) {
  data = a.data;
  subnet = a.subnet;
  env = a.env;
  actions = a.actions ? new ptrlist<analysis>(*a.actions) : nullptr;
  type = a.type;
  runs = a.runs;
  progress = a.progress;
}

void analysis::addAnalysis(analysis *a) {
  if (!actions) {
    actions = new ptrlist<analysis>();
  }
  actions->push_front(a);
}

void analysis::delAnalysis(analysis *a) {
  if (actions != nullptr) {
    actions->remove(a);
  }
}

/* Creates a sweep object depending on the analysis's properties. */
sweep *analysis::createSweep(const std::string &n) {
  sweep *swp = nullptr;
  // get type of sweep
  const char *const type = getPropertyString("Type");

  // linearly or logarithmically stepped sweeps
  if (!strcmp(type, "lin") || !strcmp(type, "log")) {
    double start = getPropertyDouble("Start");
    double stop = getPropertyDouble("Stop");
    int points = getPropertyInteger("Points");
    if (!strcmp(type, "lin")) {
      swp = new linsweep(n);
      ((linsweep *)swp)->create(start, stop, points);
    } else if (!strcmp(type, "log")) {
      swp = new logsweep(n);
      ((logsweep *)swp)->create(start, stop, points);
    }
  }

  // lists of values
  else if (!strcmp(type, "list")) {
    vector *values = getPropertyVector("Values");
    int points = values->getSize();
    swp = new lstsweep(n);
    ((lstsweep *)swp)->create(points);
    for (int i = 0; i < values->getSize(); i++)
      swp->set(i, real(values->get(i)));
  }

  // constant value
  else if (!strcmp(type, "const")) {
    double val = getPropertyDouble("Values");
    swp = new consweep(n);
    ((consweep *)swp)->create(1);
    swp->set(0, val);
  }

  swp->setParent(this);
  return swp;
}

/* Saves the given variable into the dataset.  Creates the dataset vector if necessary. */
void analysis::saveVariable(const std::string &n, nr_complex_t z, vector *f) {
  vector *d;
  if ((d = data->findVariable(n)) == nullptr) {
    d = new vector(n);
    if (f != nullptr) {
      d->setDependencies(new strlist());
      d->getDependencies()->add(f->getName());
    }
    d->setOrigin(getName());
    data->addVariable(d);
  }
  d->add(z);
}

} // namespace qucs
