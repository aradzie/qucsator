/*
 * dataset.cpp - dataset class implementation
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

#include <cstring>

#include "complex.h"
#include "dataset.h"
#include "logging.h"
#include "object.h"
#include "strlist.h"
#include "vector.h"

#include "check_citi.h"
#include "check_csv.h"
#include "check_dataset.h"
#include "check_mdl.h"
#include "check_touchstone.h"
#include "check_zvr.h"

namespace qucs {

dataset::dataset() : object() {
  variables = dependencies = nullptr;
  file = nullptr;
}

dataset::dataset(char *n) : object(n) {
  variables = dependencies = nullptr;
  file = nullptr;
}

dataset::dataset(const dataset &d) : object(d) {
  file = d.file ? strdup(d.file) : nullptr;
  vector *v;
  // copy dependency vectors
  for (v = d.dependencies; v != nullptr; v = (vector *)v->getNext()) {
    addDependency(new vector(*v));
  }
  // copy variable vectors
  for (v = variables; v != nullptr; v = (vector *)v->getNext()) {
    addVariable(new vector(*v));
  }
}

dataset::~dataset() {
  vector *n, *v;
  // delete dependency vectors
  for (v = dependencies; v != nullptr; v = n) {
    n = (vector *)v->getNext();
    delete v;
  }
  // delete variable vectors
  for (v = variables; v != nullptr; v = n) {
    n = (vector *)v->getNext();
    delete v;
  }
  free(file);
}

// Adds a dependency vector to the current dataset.
void dataset::addDependency(vector *v) {
  if (dependencies)
    dependencies->setPrev(v);
  v->setNext(dependencies);
  v->setPrev(nullptr);
  dependencies = v;
}

// Removes a dependency vector from the current dataset.
void dataset::delDependency(vector *v) {
  if (dependencies == v) {
    dependencies = (vector *)v->getNext();
    if (dependencies)
      dependencies->setPrev(nullptr);
  } else {
    vector *next = (vector *)v->getNext();
    vector *prev = (vector *)v->getPrev();
    prev->setNext(next);
    if (next)
      next->setPrev(prev);
  }
  delete v;
}

/* The function adds the given list of vectors to the dependency set
   of the current dataset. */
void dataset::addDependencies(vector *v) {
  vector *next;
  for (vector *t = v; t != nullptr; t = next) {
    next = (vector *)t->getNext();
    addDependency(t);
  }
}

// Appends a dependency vector to the current dataset.
void dataset::appendDependency(vector *v) {
  vector *e;
  if (dependencies) {
    for (e = dependencies; e->getNext(); e = (vector *)e->getNext())
      ;
    v->setPrev(e);
    e->setNext(v);
  } else {
    v->setPrev(nullptr);
    dependencies = v;
  }
  v->setNext(nullptr);
}

/* Appends the given list of vectors to the dependency
   set of the current dataset. */
void dataset::appendDependencies(vector *v) {
  vector *next;
  for (vector *t = v; t != nullptr; t = next) {
    next = (vector *)t->getNext();
    appendDependency(t);
  }
}

// Adds a variable vector to the current dataset.
void dataset::addVariable(vector *v) {
  if (variables)
    variables->setPrev(v);
  v->setNext(variables);
  v->setPrev(nullptr);
  variables = v;
}

// Removes a variable vector from the current dataset.
void dataset::delVariable(vector *v) {
  if (variables == v) {
    variables = (vector *)v->getNext();
    if (variables)
      variables->setPrev(nullptr);
  } else {
    vector *next = (vector *)v->getNext();
    vector *prev = (vector *)v->getPrev();
    prev->setNext(next);
    if (next)
      next->setPrev(prev);
  }
  delete v;
}

/* Adds the given list of vectors to the variable set of the current dataset. */
void dataset::addVariables(vector *v) {
  vector *next;
  for (vector *t = v; t != nullptr; t = next) {
    next = (vector *)t->getNext();
    addVariable(t);
  }
}

// Appends a variable vector to the current dataset.
void dataset::appendVariable(vector *v) {
  vector *e;
  if (variables) {
    for (e = variables; e->getNext(); e = (vector *)e->getNext())
      ;
    v->setPrev(e);
    e->setNext(v);
  } else {
    v->setPrev(nullptr);
    variables = v;
  }
  v->setNext(nullptr);
}

/* Appends the given list of vectors to the variable set of the current dataset. */
void dataset::appendVariables(vector *v) {
  vector *next;
  for (vector *t = v; t != nullptr; t = next) {
    next = (vector *)t->getNext();
    appendVariable(t);
  }
}

/* Applies the dependency string list of the given vector to the list of vectors
 * appended to this vector. */
void dataset::applyDependencies(vector *v) {
  strlist *deps = v->getDependencies();
  if (deps != nullptr) {
    vector *next;
    for (vector *t = (vector *)v->getNext(); t != nullptr; t = next) {
      next = (vector *)t->getNext();
      if (t->getDependencies() == nullptr) {
        t->setDependencies(new strlist(*deps));
      }
    }
  }
}

/* Returns the dataset vector (both independent and dependent) with the given origin.
 * It returns nullptr if there is no such vector. */
vector *dataset::findOrigin(char *n) {
  vector *v;
  for (v = variables; v != nullptr; v = (vector *)v->getNext()) {
    char *origin = v->getOrigin();
    if (origin != nullptr && n != nullptr && !strcmp(n, origin))
      return v;
  }
  for (v = dependencies; v != nullptr; v = (vector *)v->getNext()) {
    char *origin = v->getOrigin();
    if (origin != nullptr && n != nullptr && !strcmp(n, origin))
      return v;
  }
  return nullptr;
}

/* Assigns dependency entries to variable vectors which do have the specified origin. */
void dataset::assignDependency(const char *const origin, const char *const depvar) {
  for (vector *v = variables; v != nullptr; v = (vector *)v->getNext()) {
    char *n = v->getOrigin();
    if (n != nullptr && origin != nullptr && !strcmp(origin, n)) {
      strlist *deplist = v->getDependencies();
      if (deplist != nullptr) {
        if (!deplist->contains(depvar)) {
          deplist->append(depvar);
        }
      } else {
        deplist = new strlist();
        deplist->add(depvar);
        v->setDependencies(deplist);
      }
    }
  }
}

// Return non-zero if the given vector is an independent variable vector.
int dataset::isDependency(vector *dep) {
  for (vector *v = dependencies; v != nullptr; v = (vector *)v->getNext())
    if (v == dep)
      return 1;
  return 0;
}

// Return non-zero if the given vector is a dependent variable vector.
int dataset::isVariable(vector *var) {
  for (vector *v = variables; v != nullptr; v = (vector *)v->getNext())
    if (v == var)
      return 1;
  return 0;
}

/* Goes through the list of dependencies in the dataset
   and returns the vector specified by the given name.  Otherwise the
   function returns nullptr. */
vector *dataset::findDependency(const char *n) {
  for (vector *v = dependencies; v != nullptr; v = (vector *)v->getNext()) {
    if (!strcmp(v->getName(), n))
      return v;
  }
  return nullptr;
}

/* Goes through the list of variables in the dataset and
   returns the vector specified by the given name.  If there is no
   such variable registered the function returns nullptr. */
vector *dataset::findVariable(const std::string &name) {
  for (vector *v = variables; v != nullptr; v = (vector *)v->getNext()) {
    if (!strcmp(v->getName(), name.c_str()))
      return v;
  }
  return nullptr;
}

// Returns the number of variable vectors.
int dataset::countVariables() {
  int count = 0;
  for (vector *v = variables; v != nullptr; v = (vector *)v->getNext())
    count++;
  return count;
}

// Returns the number of dependency vectors.
int dataset::countDependencies() {
  int count = 0;
  for (vector *v = dependencies; v != nullptr; v = (vector *)v->getNext())
    count++;
  return count;
}

// Returns the current output file name.
char *dataset::getFile() { return file; }

/* Sets the current output file name.  The file name is used during
   the print functionality of the dataset class. */
void dataset::setFile(const char *f) {
  free(file);
  file = f ? strdup(f) : nullptr;
}

/* Prints the current dataset representation either to
   the specified file name (given by the function setFile()) or to
   stdout if there is no such file name given. */
void dataset::print() {
  FILE *f = stdout;

  // open file for writing
  if (file) {
    if ((f = fopen(file, "w")) == nullptr) {
      logprint(LOG_ERROR, "cannot create file `%s': %s\n", file, strerror(errno));
      return;
    }
  }

  // print header
  fprintf(f, "<Qucs Dataset>\n");

  // print dependencies
  for (vector *d = dependencies; d != nullptr; d = (vector *)d->getNext()) {
    printDependency(d, f);
  }

  // print variables
  for (vector *v = variables; v != nullptr; v = (vector *)v->getNext()) {
    if (v->getDependencies() != nullptr)
      printVariable(v, f);
    else
      printDependency(v, f);
  }

  // close file if necessary
  if (file)
    fclose(f);
}

/* Prints the given vector as independent dataset vector into the
   given file descriptor. */
void dataset::printDependency(vector *v, FILE *f) {
  // print data header
  fprintf(f, "<indep %s %d>\n", v->getName(), v->getSize());
  // print data itself
  printData(v, f);
  // print data footer
  fprintf(f, "</indep>\n");
}

/* Prints the given vector as dependent dataset vector into the given
   file descriptor. */
void dataset::printVariable(vector *v, FILE *f) {
  // print data header
  fprintf(f, "<dep %s", v->getName());
  if (v->getDependencies() != nullptr) {
    for (strlistiterator it(v->getDependencies()); *it; ++it)
      fprintf(f, " %s", *it);
  }
  fprintf(f, ">\n");

  // print data itself
  printData(v, f);

  // print data footer
  fprintf(f, "</dep>\n");
}

/* A helper routine for the print() functionality of
   the dataset class.  It prints the data items of the given vector
   object to the given output stream. */
void dataset::printData(vector *v, FILE *f) {
  for (int i = 0; i < v->getSize(); i++) {
    nr_complex_t c = v->get(i);
    if (imag(c) == 0.0) {
      fprintf(f,
              "  %+."
              "20"
              "e\n",
              real(c));
    } else {
      fprintf(f,
              "  %+."
              "20"
              "e%cj%."
              "20"
              "e\n",
              real(c), imag(c) >= 0.0 ? '+' : '-', fabs(imag(c)));
    }
  }
}

/* Reads a full dataset from the given file and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  dataset_in = f;
  dataset_restart(dataset_in);
  if (dataset_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (dataset_result != nullptr) {
    if (dataset_check(dataset_result) != 0) {
      fclose(f);
      delete dataset_result;
      return nullptr;
    }
  }
  fclose(f);
  dataset_lex_destroy();
  dataset_result->setFile(file);
  return dataset_result;
}

/* Reads a full dataset from the given touchstone file and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load_touchstone(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  touchstone_in = f;
  touchstone_restart(touchstone_in);
  if (touchstone_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (touchstone_check() != 0) {
    fclose(f);
    return nullptr;
  }
  fclose(f);
  touchstone_lex_destroy();
  touchstone_result->setFile(file);
  return touchstone_result;
}

/* Reads a full dataset from the given CSV file and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load_csv(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  csv_in = f;
  csv_restart(csv_in);
  if (csv_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (csv_check() != 0) {
    fclose(f);
    return nullptr;
  }
  fclose(f);
  csv_lex_destroy();
  csv_result->setFile(file);
  return csv_result;
}

/* Reads a full dataset from the given CITIfile and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load_citi(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  citi_in = f;
  citi_restart(citi_in);
  if (citi_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (citi_check() != 0) {
    fclose(f);
    return nullptr;
  }
  fclose(f);
  citi_lex_destroy();
  citi_result->setFile(file);
  return citi_result;
}

/* Reads a full dataset from the given ZVR file and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load_zvr(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  zvr_in = f;
  zvr_restart(zvr_in);
  if (zvr_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (zvr_check() != 0) {
    fclose(f);
    return nullptr;
  }
  fclose(f);
  zvr_lex_destroy();
  if (zvr_result)
    zvr_result->setFile(file);
  return zvr_result;
}

/* Read a full dataset from the given MDL file and returns it.
 * On failure the function emits appropriate error messages and returns nullptr. */
dataset *dataset::load_mdl(const char *file) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    logprint(LOG_ERROR, "error loading `%s': %s\n", file, strerror(errno));
    return nullptr;
  }
  mdl_in = f;
  mdl_restart(mdl_in);
  if (mdl_parse() != 0) {
    fclose(f);
    return nullptr;
  }
  if (mdl_check() != 0) {
    fclose(f);
    return nullptr;
  }
  fclose(f);
  mdl_lex_destroy();
  if (mdl_result)
    mdl_result->setFile(file);
  return mdl_result;
}

} // namespace qucs
