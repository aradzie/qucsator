/*
 * property.cpp - generic property class implementation
 *
 * Copyright (C) 2003-2009 Stefan Jahn <stefan@lkcc.org>
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

#include "property.h"
#include "complex.h"
#include "variable.h"

namespace qucs {

using namespace eqn;

property::property() : str() {
  type = PROPERTY_UNKNOWN;
  value = 0.0;
  var = nullptr;
  def = false;
}

property::~property() {
#if 0 /* FIXME: do this at another code location */
  if (type == PROPERTY_VAR) {
    constant * c = var->getConstant ();
    if (c->getType () == TAG_VECTOR) {
      delete c;
      delete var;
    }
  }
#endif
}

// Short macro in order to obtain the correct constant value.
#define D(con) ((constant *)(con))->d
#define S(con) ((constant *)(con))->s
#define V(con) ((constant *)(con))->v

// Returns the property value as vector.
qucs::vector *property::getVector() const {
  if (var != nullptr) {
    if (var->getType() == VAR_CONSTANT)
      return V(var->getConstant());
    if (var->getType() == VAR_REFERENCE)
      return V(var->getReference()->getResult());
  }
  return nullptr;
}

// Returns the property value as string.
const char *property::getString() const {
  if (var != nullptr)
    return S(var->getConstant());
  return str.c_str();
}

// Returns the property's reference if it is a variable.
const char *property::getReference() const {
  if (var != nullptr)
    return var->getName();
  return str.c_str();
}

// Returns the property value as double.
double property::getDouble() const {
  if (var != nullptr) {
    if (var->getType() == VAR_CONSTANT)
      return D(var->getConstant());
    else if (var->getType() == VAR_REFERENCE)
      return D(var->getReference()->getResult());
  }
  return value;
}

// Returns the property value as integer.
int property::getInteger() const {
  if (var != nullptr)
    return (int)std::floor(D(var->getConstant()));
  return (int)std::floor(value);
}

// Sets the property value being a double.
void property::set(const double val) {
  type = PROPERTY_DOUBLE;
  value = val;
}

// Sets the property value being an integer.
void property::set(const int val) {
  type = PROPERTY_INT;
  value = val;
}

// Sets the property value being a variable.
void property::set(variable *const val) {
  type = PROPERTY_VAR;
  var = val;
}

// Sets the property value being a string.
void property::set(const std::string &val) {
  type = PROPERTY_STR;
  this->str = val;
}

std::string property::toString() const {
  switch (type) {
  case PROPERTY_UNKNOWN:
    return "(no such type)";
  case PROPERTY_INT:
    return std::to_string(std::floor(value));
  case PROPERTY_STR:
    return std::string(this->str);
  case PROPERTY_DOUBLE:
    return std::to_string(value);
  case PROPERTY_VAR:
    return var->getName();
  }
  return "";
}

} // namespace qucs
