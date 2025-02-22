/*
 * logging.c - logging facility class implementation
 *
 * Copyright (C) 2003, 2004, 2005 Stefan Jahn <stefan@lkcc.org>
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"

/* Both of the log level dependent FILE streams. */
FILE *file_status = NULL;
FILE *file_error = NULL;

static int indent_level = 0;

void log_indent() { indent_level = indent_level + 1; }

void log_dedent() { indent_level = indent_level > 0 ? indent_level - 1 : 0; }

void logprint(const int level, const char *format, ...) {
  FILE *f = level == LOG_STATUS ? file_status : file_error;
  if (f != NULL) {
    for (int i = 0; i < indent_level; i++) {
      fputs("  ", f);
    }
    va_list args;
    va_start(args, format);
    vfprintf(f, format, args);
    va_end(args);
    fflush(f);
  }
}

/* Initialization of the logging interface. */
void loginit(void) { file_error = file_status = stderr; }

/* Customize logging. */
void redirect_status_to_stdout() { file_status = stdout; }
