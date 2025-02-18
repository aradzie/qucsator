/*
 * receiver.h - EMI receiver class definitions
 *
 * Copyright (C) 2009 Dirk Schaefer <schad@5pm.de>
 * Copyright (C) 2009 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __RECEIVER_H__
#define __RECEIVER_H__

#include <cstdint>

namespace qucs {

namespace emi {

  // receiver setting
  struct settings {
    double start;
    double stop;
    double stepsize;
    double bandwidth;
  };

  // internal helper functions
  int32_t nearestbin32 (int);
  double f_2ndorder (double, double, double);
  double f_gauss (double, double, double);
  double f_ideal (double, double, double);
  qucs::vector * receiver (double *, double, int);

  // external functionality
  qucs::vector * receiver (qucs::vector *, qucs::vector *, int len = -1);

} // namespace emi

} // namespace qucs

#endif /* __RECEIVER_H__ */
