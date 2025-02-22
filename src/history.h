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
#include <vector>
#include <utility>

namespace qucs {

class history
{
public:
  /*! default constructor */
  history ():
    sign(false),
    age(0),
    values(std::make_shared<std::vector<double>>()),
    t(std::make_shared<std::vector<double>>())
  {};

  /*! The copy constructor creates a new instance based on the given
      history object. */
  history (const history &h)
  {
      this->age = h.age;
      this->t = std::make_shared<std::vector<double>>(*(h.t));
      this->values = std::make_shared<std::vector<double>>(*(h.values));
  }

  /*! The function appends the given value to the history. */
  void push_back (const double val) {
    this->values->push_back(val);
    if (this->values != this->t)
      this->drop ();
  }

  /* This function drops the most recent n values in the history. */
  void truncate (const std::size_t n) = delete;
  void resize(const std::size_t n)
  {
    this->t->resize (n);
    this->values->resize (n);
  }

  std::size_t size (void) const
  {
    return t->size ();
  }

  void setAge (const double a) { this->age = a; }
  double getAge (void) const { return this->age; }

  // apply history
  void apply (const history & h) {
    this->t = h.t;
  }

  //! Returns the last (youngest) time value in the history
  double last (void) const {
    return this->t->empty() ? 0.0 : this->t->back();
  }

  //! Returns the first (oldest) time value in the history.
  double first (void) const {
    return this->t->empty() ? 0.0 : (*this->t)[leftidx ()];
  }

  // Returns left-most valid index into the time value vector.
  unsigned int leftidx (void) const {
    int ts = this->t->size ();
    int vs = this->values->size ();
    return ts - vs > 0 ? ts - vs : 0;
  }

  /*! Returns number of unused values (time value vector shorter than
   value vector). */
  std::size_t unused (void) {
    int ts = t->size ();
    int vs = values->size ();
    return vs - ts > 0 ? vs - ts : 0;
  }

  //! Returns the duration of the history.
  double duration(void) const {
     return last () - first ();
  }

  void truncate (const double);

  void drop (void);
  void self (void) { this->t = this->values; }

  double interpol (double, int, bool);
  double nearest (double, bool interpolate = true);
  int seek (double, int, int, double&, int);

  double getTfromidx (const int idx)  {
    return this->t == NULL ? 0.0 : (*this->t)[idx];
  }
  double getValfromidx (const int idx) {
    return this->values == NULL ? 0.0 : (*this->values)[idx];
  }

 private:
  bool sign;
  double age;
  std::shared_ptr<std::vector<double>> values;
  std::shared_ptr<std::vector<double>> t;
};

} // namespace qucs

#endif /* __HISTORY_H__ */
