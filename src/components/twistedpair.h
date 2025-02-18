/*
 * twistedpair.h - twisted pair line class definitions
 *
 * Copyright (C) 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __TWISTEDPAIR_H__
#define __TWISTEDPAIR_H__

class twistedpair : public qucs::circuit
{
 public:
  CREATOR (twistedpair);
  void initSP (void);
  void calcSP (double);
  void calcNoiseSP (double);
  void initDC (void);
  void initAC (void);
  void calcAC (double);
  void calcNoiseAC (double);
  void initTR (void);
  void saveCharacteristics (double);

 private:
  void calcPropagation (double);
  double calcLoss (double);
  double calcLength (void);

 private:
  double zl, ereff, alpha, beta, len, angle;
};

#endif /* __TWISTEDPAIR_H__ */
