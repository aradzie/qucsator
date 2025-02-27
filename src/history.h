/*
 * history.h - history class definitions
 *
 * Copyright (C) 2006, 2007 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __HISTORY_H__
#define __HISTORY_H__

#include <memory>
#include <utility>
#include <vector>

namespace qucs {

class history {
public:
  history()
      : sign(false), age(0), values(std::make_shared<std::vector<double>>()),
        t(std::make_shared<std::vector<double>>()) {};

  history(const history &h) = delete;

  /* Appends the given value to the history. */
  void push_back(const double val) {
    this->values->push_back(val);
    if (this->values != this->t) {
      this->drop();
    }
  }

  /* Drops the most recent n values in the history. */
  void truncate(const std::size_t n) = delete;
  void resize(const std::size_t n) {
    this->t->resize(n);
    this->values->resize(n);
  }

  std::size_t size() const { return t->size(); }

  void setAge(const double a) { this->age = a; }
  double getAge() const { return this->age; }

  void apply(const history &h) { this->t = h.t; }

  // Returns the last (youngest) time value in the history
  double last() const { return this->t->empty() ? 0.0 : this->t->back(); }

  // Returns the first (oldest) time value in the history.
  double first() const { return this->t->empty() ? 0.0 : (*this->t)[leftidx()]; }

  // Returns left-most valid index into the time value vector.
  unsigned int leftidx() const {
    int ts = this->t->size();
    int vs = this->values->size();
    return ts - vs > 0 ? ts - vs : 0;
  }

  /* Returns number of unused values (time value vector shorter than value vector). */
  std::size_t unused() {
    int ts = t->size();
    int vs = values->size();
    return vs - ts > 0 ? vs - ts : 0;
  }

  // Returns the duration of the history.
  double duration() const { return last() - first(); }

  void truncate(const double);

  void drop();
  void self() { this->t = this->values; }

  double interpol(double, int, bool);
  double nearest(double, bool interpolate = true);
  int seek(double, int, int, double &, int);

  double getTfromidx(const int idx) { return this->t == nullptr ? 0.0 : (*this->t)[idx]; }
  double getValfromidx(const int idx) {
    return this->values == nullptr ? 0.0 : (*this->values)[idx];
  }

private:
  bool sign;
  double age;
  std::shared_ptr<std::vector<double>> values;
  std::shared_ptr<std::vector<double>> t;
};

} // namespace qucs

#endif /* __HISTORY_H__ */
