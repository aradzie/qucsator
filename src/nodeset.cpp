/*
 * nodeset.cpp - node set class implementation
 *
 * Copyright (C) 2004, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "netdefs.h"
#include "nodeset.h"
#include "object.h"

namespace qucs {

nodeset::nodeset() {
  name = nullptr;
  value = 0.0;
  next = nullptr;
}

nodeset::nodeset(char *n) {
  name = n ? strdup(n) : nullptr;
  value = 0.0;
  next = nullptr;
}

nodeset::nodeset(char *n, double val) {
  name = n ? strdup(n) : nullptr;
  value = val;
  next = nullptr;
}

nodeset::~nodeset() { free(name); }

// Sets the name of the node set.
void nodeset::setName(char *n) {
  free(name);
  name = n ? strdup(n) : nullptr;
}

// Returns the name of the node set.
char *nodeset::getName(void) { return name; }

/* Goes through the chained list of the node sets and looks for a node
   set matching the given key and returns its value if possible.  If
   there is no such node set the function returns nullptr. */
nodeset *nodeset::findNodeset(char *n) {
  for (nodeset *p = this; p != nullptr; p = p->getNext()) {
    if (!strcmp(p->getName(), n))
      return p;
  }
  return nullptr;
}

PROP_REQ[] = {
    {"U", PROP_REAL, {0, PROP_NO_STR}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
PROP_OPT[] = {PROP_NO_PROP};
struct define_t nodeset::miscdef = {
    "NodeSet", 1, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF,
};

} // namespace qucs
