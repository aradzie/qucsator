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

namespace qucs {

class pair
{
 public:
  pair () :
    p(std::string(),0.0)
  {};

  pair (const char * const s) :
  p(s != nullptr ? std::string(s) : std::string(),0.0)
  {};

  pair (const char * const s, double v) :
    p(s != nullptr ? std::string(s) : std::string(),v)
  {} ;

  pair(const std::string s, double v) : p(s,v) {};

  void setName (const char * const s) {
    p.first = s != nullptr ? std::string(s) : std::string();
  };

  const char * getName (void) const {
    return p.first.c_str();
  };

  double getValue (void) const {
    return p.second;
  }

  void setValue (const double val) {
    p.second = val;
  }

 private:
  std::pair<std::string,double> p;
};

} // namespace qucs

#endif /* __PAIR_H__ */
