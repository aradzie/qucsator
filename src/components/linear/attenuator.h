/*
 * attenuator.h - attenuator class definitions
 *
 * Copyright (C) 2003, 2004, 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __ATTENUATOR_H__
#define __ATTENUATOR_H__

class attenuator final : public qucs::circuit {
public:
  CREATOR(attenuator);
  void initSP() override;
  void calcNoiseSP(double) override;
  void initDC() override;
  void initAC() override;
  void calcNoiseAC(double) override;
  void initTR() override;
};

#endif /* __ATTENUATOR_H__ */
