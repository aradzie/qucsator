/*
 * component_id.h - global component identifier header file
 *
 * Copyright (C) 2003-2011 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __COMPONENT_ID_H__
#define __COMPONENT_ID_H__

namespace qucs {

/* Enumerate component type identifiers. */
enum circuit_type {
  CIR_UNKNOWN = -1,

  // linear helper components
  CIR_GROUND,
  CIR_OPEN,
  CIR_SHORT,
  CIR_TEE,
  CIR_CROSS,
  CIR_ITRAFO,

  // linear components
  CIR_AMPLIFIER,
  CIR_ATTENUATOR,
  CIR_BIASTEE,
  CIR_CAPACITOR,
  CIR_CAPQ,
  CIR_CIRCULAR,
  CIR_CIRCULATOR,
  CIR_COAXLINE,
  CIR_COUPLER,
  CIR_CTLINE,
  CIR_DCBLOCK,
  CIR_DCFEED,
  CIR_GYRATOR,
  CIR_HYBRID,
  CIR_INDQ,
  CIR_INDUCTOR,
  CIR_ISOLATOR,
  CIR_MUTUAL,
  CIR_MUTUAL2,
  CIR_MUTUALX,
  CIR_PAC,
  CIR_PHASESHIFTER,
  CIR_RECTANGULAR,
  CIR_RELAIS,
  CIR_RESISTOR,
  CIR_RFEDD,
  CIR_RLCG,
  CIR_SPDFILE,
  CIR_SPFILE,
  CIR_STRAFO,
  CIR_TAPEREDLINE,
  CIR_TLINE,
  CIR_TLINE4P,
  CIR_TRAFO,
  CIR_TSWITCH,
  CIR_TWISTEDPAIR,

  // probes
  CIR_IPROBE,
  CIR_VPROBE,
  CIR_WPROBE,

  // sources
  CIR_CCCS,
  CIR_CCVS,
  CIR_IAC,
  CIR_IDC,
  CIR_IEXP,
  CIR_IFILE,
  CIR_IPULSE,
  CIR_IRECT,
  CIR_VAC,
  CIR_VAM,
  CIR_VCCS,
  CIR_VCVS,
  CIR_VDC,
  CIR_VEXP,
  CIR_VFILE,
  CIR_VPM,
  CIR_VPULSE,
  CIR_VRECT,

  // noise sources
  CIR_VNOISE,
  CIR_INOISE,
  CIR_IINOISE,
  CIR_IVNOISE,
  CIR_VVNOISE,

  // microstrip components
  CIR_MSLINE,
  CIR_MSCORNER,
  CIR_MSMBEND,
  CIR_MSSTEP,
  CIR_MSOPEN,
  CIR_MSGAP,
  CIR_MSCOUPLED,
  CIR_MSLANGE,
  CIR_MSTEE,
  CIR_MSCROSS,
  CIR_MSVIA,
  CIR_MSRSTUB,
  CIR_BONDWIRE,
  CIR_SPIRALIND,
  CIR_CIRCULARLOOP,

  // coplanar components
  CIR_CPWLINE,
  CIR_CPWOPEN,
  CIR_CPWSHORT,
  CIR_CPWGAP,
  CIR_CPWSTEP,

  // non-linear components
  CIR_OPAMP,
  CIR_DIODE,
  CIR_JFET,
  CIR_BJT,
  CIR_MOSFET,
  CIR_EQNDEFINED,
  CIR_DIAC,
  CIR_TRIAC,
  CIR_THYRISTOR,
  CIR_TUNNELDIODE,

  // digital components
  CIR_INVERTER,
  CIR_NOR,
  CIR_OR,
  CIR_NAND,
  CIR_AND,
  CIR_XNOR,
  CIR_XOR,
  CIR_DIGISOURCE,
  CIR_BUFFER,

  // external interface components
  CIR_ECVS,

};

} // namespace qucs

#endif /* __COMPONENT_ID_H__ */
