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

#include "config.h"


#include <stdlib.h>
#include <stdarg.h>

#include "logging.h"

/* Both of the log level dependent FILE streams. */
FILE * file_status = NULL;
FILE * file_error = NULL;

/* This function prints the given messages format and the appropriate
   arguments to a FILE stream depending on the given log level. */
void logprint (int level, const char * format, ...) {
  FILE * f;
  va_list args;

  f = level == LOG_STATUS ? file_status : file_error;
  if (f != NULL) {
    va_start (args, format);
    vfprintf (f, format, args);
    va_end (args);
    fflush (f);
  }
}

/* Initialization of the logging interface. */
void loginit (void) {
  file_error = file_status = stderr;
}

/* Customize logging. */
void redirect_status_to_stdout() {
  file_status = stdout;
}
