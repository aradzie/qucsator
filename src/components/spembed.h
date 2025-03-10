/*
 * spembed.h - S-parameters embedding component class definitions
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2009 Stefan Jahn <stefan@lkcc.org>
 * Copyright (C) 2017 Qucs Team
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
 *
 */

#ifndef SPEMBED_H
#define SPEMBED_H

#include "spfile.h"

/* S-parameters embedding component takes data from an SnP file */
class spembed : public spfile, public qucs::circuit
{
 public:
  CREATOR(spembed);

  void initSP (void);
  void calcSP (double);
  void initDC (void);
  void initTR (void);
  void initAC (void);
  void calcAC (double);

  void calcNoiseSP (double);
  void calcNoiseAC (double);

  qucs::matrix correlationMatrix (double, nr_complex_t, double, qucs::matrix);
  double noiseFigure (qucs::matrix, qucs::matrix, double&, nr_complex_t&, double&);
  qucs::matrix expandNoiseMatrix (qucs::matrix, qucs::matrix);
  qucs::matrix shrinkNoiseMatrix (qucs::matrix, qucs::matrix);
  qucs::matrix calcMatrixCs (double);
};

#endif /* SPEMBED_H */
