/*
 * rfedd.h - equation defined RF device class definitions
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

#ifndef __RFEDD_H__
#define __RFEDD_H__

class rfedd : public qucs::circuit
{
 public:
  CREATOR (rfedd);
  ~rfedd ();
  void initDC (void);
  void calcDC (void);
  void initAC (void);
  void calcAC (double);
  void initSP (void);
  void calcSP (double);
  void initTR (void);
  void calcTR (double);

 private:
  void initModel (void);
  char * createVariable (const char *, int, int, bool prefix = true);
  char * createVariable (const char *, bool prefix = true);
  void setResult (void *, double);
  void setResult (void *, nr_complex_t);
  nr_complex_t getResult (void *);
  qucs::matrix calcMatrix (double);
  void updateLocals (double);
  void prepareModel (void);
  void initMNA (void);
  void calcMNA (double);

 private:
  void ** peqn;
  void * seqn;
  void * feqn;
};

#endif /* __RFEDD_H__ */
