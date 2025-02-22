/*
 * exceptionstack.cpp - exception stack class implementation
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
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

#include "exceptionstack.h"
#include "exception.h"
#include "logging.h"

using namespace qucs;

// The global exception stack.
exceptionstack qucs::estack;

exceptionstack::exceptionstack() : root(nullptr) {}

exceptionstack::~exceptionstack() {
  while (root) {
    exception *next = root->getNext();
    delete root;
    root = next;
  }
}

// Pushes a new exception onto the exception stack.
void exceptionstack::push(exception *e) {
  e->setNext(root);
  root = e;
}

/* Removes the top exception from the exception stack and returns the new top exception. */
exception *exceptionstack::pop() {
  if (root != nullptr) {
    exception *next = root->getNext();
    delete root;
    root = next;
  }
  return root;
}

// Returns the top exception.
exception *exceptionstack::top() { return root; }

/* Prints the complete exception stack and removes each exception from the stack. */
void exceptionstack::print(const char *prefix) {
  if (root) {
    logprint(LOG_ERROR, "%s%sexception stack\n", prefix ? prefix : "", prefix ? " " : "");
  }
  exception *next;
  while ((next = top()) != nullptr) {
    logprint(LOG_ERROR, "  %03d: %s\n", next->getCode(), next->getText());
    pop();
  }
}
