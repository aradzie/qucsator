/*
 * environment.h - variable environment class definitions
 *
 * Copyright (C) 2004, 2006, 2007, 2009 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__

#include <list>
#include <string>

#include "equation.h"

namespace qucs {

class variable;
class checker;
class solver;
class dataset;

/* The environment class holds information and pointers to the classes and methods
 * used to evaluate a netlist. */
class environment {
public:
  environment();
  environment(const std::string &p_name);
  environment(const environment &);
  virtual ~environment();
  void copy(const environment &);
  void setName(char *) = delete;
  void print(bool all = false) const;
  void setDefinitions(definition_t *const d) { defs = d; }
  definition_t *getDefinitions() const { return defs; }

  // variable specific functionality
  void copyVariables(variable *);
  void deleteVariables();
  void addVariable(variable *, bool pass = true);
  variable *getVariable(const char *) const;

  // equation specific functionality
  void setChecker(eqn::checker *c) { checkee = c; }
  eqn::checker *getChecker() { return checkee; }
  void setSolver(eqn::solver *s) { solvee = s; }
  eqn::solver *getSolver() { return solvee; }
  int equationChecker(int noundefined = 1) const;
  int equationSolver(dataset *);
  int runSolver();
  void equationSolver();

  // subcircuit specific
  qucs::vector getVector(const char *) const;
  void setDoubleConstant(const char *, double);
  double getDoubleConstant(const char *) const;
  void setDouble(const char *, double);
  double getDouble(const char *) const;
  void setDoubleReference(const char *, char *);
  char *getDoubleReference(const char *) const;
  void updateReferences(environment *);
  void passConstants();
  void fetchConstants();
  variable *findValue(char *);
  void setValue(char *, eqn::constant *);
  void saveResults();

  /* Adds a child to the environment. */
  inline void push_front_Child(environment *child) { children.push_front(child); }

  /* Removes a child from the environment. */
  void remove_Child(environment *child) { children.remove(child); }

  /* set the name */
  void setName(const std::string &p_name) { this->name = p_name; }

  /* Returns the name of the environment. */
  const std::string &getName() const { return this->name; }

private:
  std::string name;
  variable *root;
  eqn::checker *checkee;
  eqn::solver *solvee;
  std::list<environment *> children;
  bool iscopy;
  struct definition_t *defs;
};

} // namespace qucs

#endif /* __ENVIRONMENT_H__ */
