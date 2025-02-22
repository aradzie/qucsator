/*
 * exception.cpp - exception class implementation
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

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "exception.h"

using namespace qucs;

exception::exception() : code(EXCEPTION_UNKNOWN), data(0), text(nullptr), next(nullptr) {}

exception::exception(const int code) : code(code), data(0), text(nullptr), next(nullptr) {}

exception::~exception() { free(text); }

void exception::setText(const char *format, ...) {
  va_list args;

  const auto str = (char *)malloc(1024); // this should be enough
  va_start(args, format);
  vsprintf(str, format, args);
  va_end(args);

  free(text);
  text = strdup(str);
  free(str);
}
