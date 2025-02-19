/*
 * variable.h - generic variable class definitions
 *
 * Copyright (C) 2004, 2007 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#include <string>

#include "analysis.h"
#include "components/microstrip/substrate.h"
#include "equation.h"

namespace qucs {

enum variably_type {
  VAR_UNKNOWN = -1,
  VAR_CONSTANT,
  VAR_REFERENCE,
  VAR_SUBSTRATE,
  VAR_VALUE,
  VAR_ANALYSIS
};

class substrate;
class analysis;

namespace eqn {
class equation;
class constant;
} // namespace eqn

class variable {
public:
  variable();
  variable(const char *n);
  variable(const variable &);
  virtual ~variable() = default;
  void setName(const char *const n) { name = n ? std::string(n) : std::string(); };
  const char *getName() const { return this->name.c_str(); };
  void setNext(variable *const v) { next = v; }
  variable *getNext() const { return next; }
  void setType(const int t) { type = t; }
  int getType() const { return type; }
  void setConstant(eqn::constant *const c) {
    type = VAR_CONSTANT;
    value.c = c;
  }
  eqn::constant *getConstant() const { return value.c; }
  void setReference(eqn::reference *const r) {
    type = VAR_REFERENCE;
    value.r = r;
  }
  eqn::reference *getReference() const { return value.r; }
  void setSubstrate(substrate *const s) {
    type = VAR_SUBSTRATE;
    value.s = s;
  }
  substrate *getSubstrate() { return value.s; }
  void setValue(eqn::constant *const v) {
    type = VAR_VALUE;
    value.v = v;
  }
  eqn::constant *getValue() { return value.v; }
  void setAnalysis(analysis *const a) {
    type = VAR_ANALYSIS;
    value.a = a;
  }
  analysis *getAnalysis() const { return this->value.a; }
  const char *toString();
  void setPassing(const bool p) { this->pass = p; }
  bool getPassing() const { return this->pass; }

private:
  std::string name;
  bool pass;
  int type;
  union value_t {
    eqn::constant *c;
    eqn::reference *r;
    substrate *s;
    eqn::constant *v;
    analysis *a;
  } value;
  variable *next;
};

} // namespace qucs

#endif /* __VARIABLE_H__ */
