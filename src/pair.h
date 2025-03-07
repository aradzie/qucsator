/*
 * pair.h - key/value pair class definitions
 *
 * Copyright (C) 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __PAIR_H__
#define __PAIR_H__

#include <string>

namespace qucs {

class pair {
public:
  pair(const std::string s, double v) : p(s, v) {};
  const char *getName() const { return p.first.c_str(); };
  double getValue() const { return p.second; }
  void setValue(const double value) { p.second = value; }

private:
  std::pair<std::string, double> p;
};

} // namespace qucs

#endif /* __PAIR_H__ */
