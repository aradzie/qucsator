/*
 * object.cpp - generic object class implementation
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include "object.h"
#include "property.h"
#include "variable.h"

namespace qucs {

/* Adds a property consisting of a key and a string value to the object. */
void object::addProperty(const std::string &n, const char *const val, const bool def) {
  property p;
  p.set(val);
  p.setDefault(def);
  props.insert({{n, p}});
}

/* Sets the specified property consisting of a key and a string value in the object. */
void object::setProperty(const std::string &n, const char *const val) {
  auto it = props.find(n);
  if (it != props.end()) {
    (*it).second.set(val);
  } else {
    addProperty(n, val);
  }
}

/* Adds a property consisting of a key and a double value to the object. */
void object::addProperty(const std::string &n, const double val, const bool def) {
  property p;
  p.set(val);
  p.setDefault(def);
  props.insert({{n, p}});
}

/* Sets the specified property consisting of a key and a double value in the object. */
void object::setProperty(const std::string &n, const double val) {
  auto it = props.find(n);
  if (it != props.end()) {
    (*it).second.set(val);
  } else {
    addProperty(n, val);
  }
}

/* Sets the specified property consisting of a key and a double value in the object.
 * The property is marked a scalability property. */
void object::setScaledProperty(const std::string &n, const double val) {
  setProperty("Scaled:" + n, val);
}

/* Adds a property consisting of a key and a variable value to the object. */
void object::addProperty(const std::string &n, variable *const val, const bool def) {
  property p;
  p.set(val);
  p.setDefault(def);
  props.insert({{n, p}});
}

/* Returns the requested property value which has been previously added as its vector
 * representation. If there is no such property the function returns nullptr. */
qucs::vector *object::getPropertyVector(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return (*it).second.getVector();
  }
  return nullptr;
}

/* Returns the requested property value which has been previously added as its text representation.
   If there is no such property the function returns nullptr. */
const char *object::getPropertyString(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return (*it).second.getString();
  }
  return nullptr;
}

/* Returns the requested property reference variable name.
 * If there is no such property the function returns nullptr. */
const char *object::getPropertyReference(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return (*it).second.getReference();
  }
  return nullptr;
}

/* Returns the requested property value which has been previously added as its double
 * representation. If there is no such property the function returns zero. */
double object::getPropertyDouble(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return (*it).second.getDouble();
  }
  return 0.0;
}

/* Returns the requested (scalability) property value which has been previously added.
 * If there is no such scaled property the function returns the standard property or zero. */
double object::getScaledProperty(const std::string &n) const {
  std::string txt = "Scaled:" + n;
  const auto &it = props.find(txt);
  if (it != props.end()) {
    return (*it).second.getDouble();
  }
  return this->getPropertyDouble(n);
}

/* Returns the requested property value which has been previously added as its integer
 * representation. If there is no such property the function returns zero. */
int object::getPropertyInteger(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return (*it).second.getInteger();
  }
  return 0;
}

/* Checks whether the object has got a certain property value.
 * If so it returns non-zero, otherwise it returns zero. */
bool object::hasProperty(const std::string &n) const { return props.find(n) != props.end(); }

/* Checks whether the object has got a certain property value and if this has its default value.
 * If so it returns  non-zero, otherwise it returns zero. */
bool object::isPropertyGiven(const std::string &n) const {
  const auto &it = props.find(n);
  if (it != props.end()) {
    return !(*it).second.isDefault();
  }
  return false;
}

// Returns the number of properties in the object.
int object::countProperties() const { return props.size(); }

// Returns a text representation of the objects properties.
const char *object::propertyList() const {
  std::string ptxt;
  for (auto it = props.cbegin(); it != props.cend(); ++it) {
    std::string n = it->first;
    std::string val = it->second.toString();
    std::string text = n + "=\"" + val + "\"";
    ptxt += text;
  }
  return ptxt.c_str();
}

} // namespace qucs
