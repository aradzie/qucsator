/*
 * tvector.cpp - simple vector template class implementation
 *
 * Copyright (C) 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include <cassert>
#include <vector>

#include "complex.h"
#include "tvector.h"

namespace qucs {

// Returns the tvector element at the given position.
template <class nr_type_t> nr_type_t tvector<nr_type_t>::get(int i) { return data.at(i); }

// Sets the tvector element at the given position.
template <class nr_type_t> void tvector<nr_type_t>::set(int i, nr_type_t z) { data.at(i) = z; }

// Sets all the tvector elements to the given value.
template <class nr_type_t> void tvector<nr_type_t>::set(nr_type_t z) {
  for (std::size_t i = 0; i < data.size(); i++) {
    data[i] = z;
  }
}

// Sets the specified tvector elements to the given value.
template <class nr_type_t> void tvector<nr_type_t>::set(nr_type_t z, int start, int stop) {
  for (std::size_t i = start; i < stop; i++) {
    data[i] = z;
  }
}

// Sets size to zero.  Does not reduce the capacity.
template <class nr_type_t> void tvector<nr_type_t>::clear() { data.clear(); }

/* Returns the number of entries with the given value deviating no more than the given epsilon. */
template <class nr_type_t> int tvector<nr_type_t>::contains(nr_type_t val, double eps) {
  int count = 0;
  for (int i = 0; i < data.size(); i++) {
    if (abs(data[i] - val) <= eps) {
      count++;
    }
  }
  return count;
}

// Copies the specified elements from the given tvector.
template <class nr_type_t> void tvector<nr_type_t>::set(tvector<nr_type_t> a, int start, int stop) {
  for (int i = start; i < stop; i++) {
    (*data)[i] = a.get(i);
  }
}

// Swaps the given rows with each other.
template <class nr_type_t> void tvector<nr_type_t>::exchangeRows(int r1, int r2) {
  assert(r1 >= 0 && r2 >= 0 && r1 < data.size() && r2 < data.size());
  nr_type_t s = data[r1];
  data[r1] = data[r2];
  data[r2] = s;
}

// Addition.
template <class nr_type_t>
tvector<nr_type_t> operator+(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, a.get(i) + b.get(i));
  }
  return res;
}

// Intrinsic vector addition.
template <class nr_type_t> tvector<nr_type_t> tvector<nr_type_t>::operator+=(tvector<nr_type_t> a) {
  assert(a.size() == data.size());
  std::vector<nr_type_t> *src = a.getData();
  std::vector<nr_type_t> *dst = data;
  for (int i = 0; i < data.size(); i++) {
    (*dst)[i] += (*src)[i];
  }
  return *this;
}

// Subtraction.
template <class nr_type_t>
tvector<nr_type_t> operator-(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, a.get(i) - b.get(i));
  }
  return res;
}

// Intrinsic vector subtraction.
template <class nr_type_t> tvector<nr_type_t> tvector<nr_type_t>::operator-=(tvector<nr_type_t> a) {
  assert(a.size() == data.size());
  std::vector<nr_type_t> *src = a.getData();
  std::vector<nr_type_t> *dst = data;
  for (int i = 0; i < data.size(); i++) {
    (*dst)[i] -= (*src)[i];
  }
  return *this;
}

// Intrinsic scalar multiplication.
template <class nr_type_t> tvector<nr_type_t> tvector<nr_type_t>::operator*=(double s) {
  std::vector<nr_type_t> *dst = data;
  for (int i = 0; i < data.size(); i++) {
    (*dst)[i] *= s;
  }
  return *this;
}

// Intrinsic scalar division.
template <class nr_type_t> tvector<nr_type_t> tvector<nr_type_t>::operator/=(double s) {
  std::vector<nr_type_t> *dst = data;
  for (int i = 0; i < data.size(); i++) {
    (*dst)[i] /= s;
  }
  return *this;
}

// Scalar multiplication.
template <class nr_type_t> tvector<nr_type_t> operator*(double s, tvector<nr_type_t> a) {
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, s * a.get(i));
  }
  return res;
}

template <class nr_type_t> tvector<nr_type_t> operator*(tvector<nr_type_t> a, double s) {
  return s * a;
}

// Vector multiplication (element by element).
template <class nr_type_t>
tvector<nr_type_t> operator*(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, a.get(i) * b.get(i));
  }
  return res;
}

// Computes the scalar product of two vectors.
template <class nr_type_t> nr_type_t scalar(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  nr_type_t n = 0;
  for (int i = 0; i < a.size(); i++) {
    n += a.get(i) * b.get(i);
  }
  return n;
}

// Constant assignment operation.
template <class nr_type_t> tvector<nr_type_t> tvector<nr_type_t>::operator=(const nr_type_t val) {
  for (int i = 0; i < data.size(); i++) {
    (*data)[i] = val;
  }
  return *this;
}

// Returns the sum of the vector elements.
template <class nr_type_t> nr_type_t sum(tvector<nr_type_t> a) {
  nr_type_t res = 0;
  for (int i = 0; i < a.size(); i++) {
    res += a.get(i);
  }
  return res;
}

// Vector negation.
template <class nr_type_t> tvector<nr_type_t> operator-(tvector<nr_type_t> a) {
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, -a.get(i));
  }
  return res;
}

// Vector less comparison.
template <class nr_type_t> bool operator<(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  for (int i = 0; i < a.size(); i++) {
    if (a.get(i) >= b.get(i)) {
      return false;
    }
  }
  return true;
}

// Vector greater comparison.
template <class nr_type_t> bool operator>(tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert(a.size() == b.size());
  for (int i = 0; i < a.size(); i++) {
    if (a.get(i) <= b.get(i)) {
      return false;
    }
  }
  return true;
}

// Scalar addition.
template <class nr_type_t> tvector<nr_type_t> operator+(nr_type_t s, tvector<nr_type_t> a) {
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, s + a.get(i));
  }
  return res;
}

template <class nr_type_t> tvector<nr_type_t> operator+(tvector<nr_type_t> a, nr_type_t s) {
  return s + a;
}

// Mean square norm.
template <class nr_type_t> double norm(tvector<nr_type_t> a) {
#if 0
  double k = 0;
  for (int i = 0; i < a.size (); i++) k += norm (a.get (i));
  return n;
#else
  double scale = 0, n = 1, x, ax;
  for (int i = 0; i < a.size(); i++) {
    if ((x = real(a(i))) != 0) {
      ax = fabs(x);
      if (scale < ax) {
        x = scale / ax;
        n = 1 + n * x * x;
        scale = ax;
      } else {
        x = ax / scale;
        n += x * x;
      }
    }
    if ((x = imag(a(i))) != 0) {
      ax = fabs(x);
      if (scale < ax) {
        x = scale / ax;
        n = 1 + n * x * x;
        scale = ax;
      } else {
        x = ax / scale;
        n += x * x;
      }
    }
  }
  return scale * scale * n;
#endif
}

// Maximum norm.
template <class nr_type_t> double maxnorm(tvector<nr_type_t> a) {
  double nMax = 0;
  for (int i = 0; i < a.size(); i++) {
    double n = norm(a.get(i));
    if (n > nMax) {
      nMax = n;
    }
  }
  return nMax;
}

// Conjugate vector.
template <class nr_type_t> tvector<nr_type_t> conj(tvector<nr_type_t> a) {
  tvector<nr_type_t> res((a.size()));
  for (int i = 0; i < a.size(); i++) {
    res.set(i, conj(a.get(i)));
  }
  return res;
}

// Checks validity of vector.
template <class nr_type_t> int tvector<nr_type_t>::isFinite() {
  for (int i = 0; i < data.size(); i++) {
    if (!std::isfinite(real((*data)[i]))) {
      return 0;
    }
  }
  return 1;
}

// Reorders the vector according to the given index array.
template <class nr_type_t> void tvector<nr_type_t>::reorder(int *idx) {
  tvector<nr_type_t> old = *this;
  for (int i = 0; i < data.size(); i++) {
    (*data)[i] = old.get(idx[i]);
  }
}

#ifdef DEBUG
template <class nr_type_t> void tvector<nr_type_t>::print() {
  for (int r = 0; r < data.size(); r++) {
    fprintf(stderr, "%+.2e%+.2ei\n", (double)real(get(r)), (double)imag(get(r)));
  }
}
#endif /* DEBUG */

} // namespace qucs
