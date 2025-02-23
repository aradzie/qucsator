/*
 * indq.h - Lossy inductor class definition
 *
 * Copyright (C) 2015 Andres Martinez-Mera <andresmartinezmera@gmail.com>
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

#ifndef INDQ_H
#define INDQ_H

class indq final : public qucs::circuit {
public:
  CREATOR(indq);
  void calcSP(double) override;
  void calcNoiseSP(double) override;
  void initDC() override;
  void calcDC() override;
  void initAC() override;
  void initSP() override;
  void calcAC(double) override;
  void calcNoiseAC(double) override;

private:
  void calcZs(double);
  nr_complex_t Zs;
};

#endif /* __indq_H__ */
