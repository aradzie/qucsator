/*
 * module.h - module class definitions
 *
 * Copyright (C) 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __MODULE_H__
#define __MODULE_H__

#include "hash.h"

namespace qucs {

class circuit;
class analysis;

// function typedefs for circuits and analyses
typedef circuit *(*circuit_creator_t)();
typedef analysis *(*analysis_creator_t)();
typedef struct define_t *(*circuit_definer_t)();
typedef struct define_t *(*analysis_definer_t)();
typedef struct define_t *(*misc_definer_t)();

class module {
public:
  module();
  ~module();

  static void registerModule(circuit_definer_t, circuit_creator_t);
  static void registerModules();
  static struct define_t *getModule(char *);

private:
  static void registerModule(analysis_definer_t, analysis_creator_t);
  static void registerModule(struct define_t *);
  static void registerModule(misc_definer_t);
  static void registerModule(const char *, module *);

public:
  static qucs::hash<module> modules;

public:
  struct define_t *definition;
  circuit_creator_t circreate;
  analysis_creator_t anacreate;
};

} // namespace qucs

#endif /* __MODULE_H__ */
