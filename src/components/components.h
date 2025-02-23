/*
 * components.h - global component header file
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

#ifndef __COMPONENTS_H__
#define __COMPONENTS_H__

// BUG: include all component headers.
// components should add to the kernel, not the other way around.

#include "component_id.h"

#include "circuit.h"
#include "complex.h"
#include "node.h"
#include "object.h"

#include "helpers/cross.h"
#include "helpers/ground.h"
#include "helpers/itrafo.h"
#include "helpers/open.h"
#include "helpers/short.h"
#include "helpers/tee.h"

#include "linear/amplifier.h"
#include "linear/attenuator.h"
#include "linear/biastee.h"
#include "linear/capacitor.h"
#include "linear/capq.h"
#include "linear/coupler.h"
#include "linear/dcblock.h"
#include "linear/dcfeed.h"
#include "linear/gyrator.h"
#include "linear/hybrid.h"
#include "linear/indq.h"
#include "linear/inductor.h"
#include "linear/isolator.h"
#include "linear/mutual.h"
#include "linear/mutual2.h"
#include "linear/mutualx.h"
#include "linear/phaseshifter.h"
#include "linear/relais.h"
#include "linear/resistor.h"
#include "linear/tswitch.h"

#include "noise/iinoise.h"
#include "noise/inoise.h"
#include "noise/ivnoise.h"
#include "noise/vnoise.h"
#include "noise/vvnoise.h"

#include "sources/cccs.h"
#include "sources/ccvs.h"
#include "sources/iac.h"
#include "sources/idc.h"
#include "sources/iexp.h"
#include "sources/ifile.h"
#include "sources/ipulse.h"
#include "sources/irect.h"
#include "sources/vac.h"
#include "sources/vam.h"
#include "sources/vccs.h"
#include "sources/vcvs.h"
#include "sources/vdc.h"
#include "sources/vexp.h"
#include "sources/vfile.h"
#include "sources/vpm.h"
#include "sources/vpulse.h"
#include "sources/vrect.h"

#include "probes/iprobe.h"
#include "probes/vprobe.h"
#include "probes/wprobe.h"

#include "nonlinear/bjt.h"
#include "nonlinear/diac.h"
#include "nonlinear/diode.h"
#include "nonlinear/eqndefined.h"
#include "nonlinear/jfet.h"
#include "nonlinear/mosfet.h"
#include "nonlinear/opamp.h"
#include "nonlinear/thyristor.h"
#include "nonlinear/triac.h"
#include "nonlinear/tunneldiode.h"

#include "microstrip/substrate.h"

#include "microstrip/bondwire.h"
#include "microstrip/circularloop.h"
#include "microstrip/mscorner.h"
#include "microstrip/mscoupled.h"
#include "microstrip/mscross.h"
#include "microstrip/msgap.h"
#include "microstrip/mslange.h"
#include "microstrip/msline.h"
#include "microstrip/msmbend.h"
#include "microstrip/msopen.h"
#include "microstrip/msrstub.h"
#include "microstrip/msstep.h"
#include "microstrip/mstee.h"
#include "microstrip/msvia.h"
#include "microstrip/spiralinductor.h"

#include "microstrip/cpwgap.h"
#include "microstrip/cpwline.h"
#include "microstrip/cpwopen.h"
#include "microstrip/cpwshort.h"
#include "microstrip/cpwstep.h"

#include "digital/and.h"
#include "digital/buffer.h"
#include "digital/digisource.h"
#include "digital/inverter.h"
#include "digital/nand.h"
#include "digital/nor.h"
#include "digital/or.h"
#include "digital/xnor.h"
#include "digital/xor.h"

#include "pac.h"
#include "circulator.h"
#include "trafo.h"
#include "strafo.h"
#include "tline.h"
#include "ctline.h"
#include "coaxline.h"
#include "circline.h"
#include "taperedline.h"
#include "rectline.h"
#include "twistedpair.h"
#include "tline4p.h"
#include "rlcg.h"
#include "spembed.h"
#include "spdeembed.h"
#include "rfedd.h"

#include "ecvs.h"

#endif /* __COMPONENTS_H__ */
