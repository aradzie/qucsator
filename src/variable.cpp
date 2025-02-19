/*
 * variable.cpp - generic variable class implementation
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

#include "variable.h"
#include "analysis.h"
#include "components/microstrip/substrate.h"
#include "equation.h"

namespace qucs {

variable::variable() : name() {
  next = nullptr;
  type = VAR_UNKNOWN;
  pass = true;
}

variable::variable(const char *const n) {
  name = n ? std::string(n) : std::string();
  next = nullptr;
  type = VAR_UNKNOWN;
  pass = true;
}

variable::variable(const variable &o) {
  this->name = o.name;
  type = o.type;
  next = o.next;
  pass = o.pass;
  value = o.value;
}

const char *variable::toString() {
  std::string text;
  const char *str = nullptr;
  char *val = nullptr;
  switch (type) {
  case VAR_UNKNOWN:
    text = "variable";
    break;
  case VAR_CONSTANT:
    str = value.c->toString();
    text = "constant: " + std::string(str);
    break;
  case VAR_VALUE:
    str = value.v->toString();
    text = "value: " + std::string(str);
    break;
  case VAR_REFERENCE:
    str = value.r->toString();
    val = value.r->getResult()->toString();
    text = "reference: " + std::string(str) + " = " + std::string(val);
    break;
  case VAR_SUBSTRATE:
    str = value.s->getName();
    text = "substrate: " + std::string(str);
    break;
  case VAR_ANALYSIS:
    str = value.a->getName();
    text = "analysis: " + std::string(str);
    break;
  default:
    text = "?variable?";
    break;
  }
  return text.c_str();
}

} // namespace qucs
