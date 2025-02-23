/*
 * module.cpp - module class implementation
 *
 * Copyright (C) 2008, 2009, 2010 Stefan Jahn <stefan@lkcc.org>
 * New models added Mike Brinson 2013 mbrin72043@yahoo.co.uk
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

#include <cstdio>
#include <iostream>
#include <list>
#include <map>

#include "analyses.h"
#include "components.h"
#include "logging.h"
#include "module.h"
#include "netdefs.h"
#include "nodeset.h"

namespace qucs {

// Global module hash.
qucs::hash<module> module::modules;

// Our global factories for making loaded circuit objects
// The factories are populated as dynamic modules are loaded
// factorycreate hold the loaded modules constructor function
std::map<std::string, creator_t *> factorycreate;
// factorydef hold the loaded modules definitions
std::map<std::string, defs_t *> factorydef;

// Constructor creates an instance of the module class.
module::module() {
  definition = nullptr;
  circreate = nullptr;
  anacreate = nullptr;
}

module::~module() {}

// Definitions which do not fit elsewhere.
static struct property_t props1[] = {
    PROP_NO_PROP,
};
static struct property_t props2[] = {
    {"Type", PROP_STR, {PROP_NO_VAL, "DEF1"}, PROP_NO_RANGE},
    PROP_NO_PROP,
};
struct define_t miscdef1 = {
    "Def", PROP_NODES, PROP_ACTION, PROP_NO_SUBSTRATE, PROP_LINEAR, props1, props1,
};
struct define_t miscdef2 = {
    "Sub", PROP_NODES, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, props1, props2,
};

// Registers an analysis object to the list of available modules.
void module::registerModule(analysis_definer_t define, analysis_creator_t create) {
  module *m = new module();
  m->definition = define();
  m->anacreate = create;
  modules.put((char *)define()->type, m);
}

// Registers a circuit object to the list of available modules.
void module::registerModule(circuit_definer_t define, circuit_creator_t create) {
  module *m = new module();
  m->definition = define();
  m->circreate = create;
  registerModule(define()->type, m);
}

// Registers a miscellaneous object to the list of available modules.
void module::registerModule(misc_definer_t define) {
  module *m = new module();
  m->definition = define();
  registerModule(define()->type, m);
}

/* Registers a miscellaneous structure just defined by a simple
   define_t structure to the list of available modules. */
void module::registerModule(struct define_t *define) {
  module *m = new module();
  m->definition = define;
  registerModule(define->type, m);
}

// Puts a module into the available module hash.
void module::registerModule(const char *type, module *m) {
  if (modules.get((char *)type) != nullptr) {
    logprint(LOG_ERROR, "module already registered: %s\n", type);
  } else {
    modules.put((char *)type, m);
  }
}

/* Returns the definition of a module specified by its type name if
   such is existing and otherwise nullptr. */
struct define_t *module::getModule(char *type) {
  module *m = modules.get(type);
  if (m != nullptr) {
    return m->definition;
  }
  return nullptr;
}

// Helper macros.
#define REGISTER_CIRCUIT(val) registerModule(&val::definition, &val::create)
#define REGISTER_ANALYSIS(val) registerModule(&val::definition, &val::create)
#define REGISTER_MISC(val) registerModule(&val::definition)

// Global static module registration.
void module::registerModules() {
  // miscellaneous
  registerModule(&miscdef1);
  registerModule(&miscdef2);
  REGISTER_MISC(nodeset);
  REGISTER_MISC(substrate);

  // circuit components
  REGISTER_CIRCUIT(amplifier);
  REGISTER_CIRCUIT(attenuator);
  REGISTER_CIRCUIT(biastee);
  REGISTER_CIRCUIT(bjt);
  REGISTER_CIRCUIT(bondwire);
  REGISTER_CIRCUIT(buffer);
  REGISTER_CIRCUIT(capacitor);
  REGISTER_CIRCUIT(capq);
  REGISTER_CIRCUIT(cccs);
  REGISTER_CIRCUIT(ccvs);
  REGISTER_CIRCUIT(circline);
  REGISTER_CIRCUIT(circularloop);
  REGISTER_CIRCUIT(circulator);
  REGISTER_CIRCUIT(coaxline);
  REGISTER_CIRCUIT(coupler);
  REGISTER_CIRCUIT(cpwgap);
  REGISTER_CIRCUIT(cpwline);
  REGISTER_CIRCUIT(cpwopen);
  REGISTER_CIRCUIT(cpwshort);
  REGISTER_CIRCUIT(cpwstep);
  REGISTER_CIRCUIT(ctline);
  REGISTER_CIRCUIT(dcblock);
  REGISTER_CIRCUIT(dcfeed);
  REGISTER_CIRCUIT(diac);
  REGISTER_CIRCUIT(digisource);
  REGISTER_CIRCUIT(diode);
  REGISTER_CIRCUIT(eqndefined);
  REGISTER_CIRCUIT(gyrator);
  REGISTER_CIRCUIT(hybrid);
  REGISTER_CIRCUIT(iac);
  REGISTER_CIRCUIT(idc);
  REGISTER_CIRCUIT(iexp);
  REGISTER_CIRCUIT(ifile);
  REGISTER_CIRCUIT(iinoise);
  REGISTER_CIRCUIT(indq);
  REGISTER_CIRCUIT(inductor);
  REGISTER_CIRCUIT(inoise);
  REGISTER_CIRCUIT(inverter);
  REGISTER_CIRCUIT(iprobe);
  REGISTER_CIRCUIT(ipulse);
  REGISTER_CIRCUIT(irect);
  REGISTER_CIRCUIT(isolator);
  REGISTER_CIRCUIT(ivnoise);
  REGISTER_CIRCUIT(jfet);
  REGISTER_CIRCUIT(logicand);
  REGISTER_CIRCUIT(logicnand);
  REGISTER_CIRCUIT(logicnor);
  REGISTER_CIRCUIT(logicor);
  REGISTER_CIRCUIT(logicxnor);
  REGISTER_CIRCUIT(logicxor);
  REGISTER_CIRCUIT(mosfet);
  REGISTER_CIRCUIT(mscorner);
  REGISTER_CIRCUIT(mscoupled);
  REGISTER_CIRCUIT(mscross);
  REGISTER_CIRCUIT(msgap);
  REGISTER_CIRCUIT(mslange);
  REGISTER_CIRCUIT(msline);
  REGISTER_CIRCUIT(msmbend);
  REGISTER_CIRCUIT(msopen);
  REGISTER_CIRCUIT(msrstub);
  REGISTER_CIRCUIT(msstep);
  REGISTER_CIRCUIT(mstee);
  REGISTER_CIRCUIT(msvia);
  REGISTER_CIRCUIT(mutual);
  REGISTER_CIRCUIT(mutual2);
  REGISTER_CIRCUIT(mutualx);
  REGISTER_CIRCUIT(opamp);
  REGISTER_CIRCUIT(pac);
  REGISTER_CIRCUIT(phaseshifter);
  REGISTER_CIRCUIT(rectline);
  REGISTER_CIRCUIT(relais);
  REGISTER_CIRCUIT(resistor);
  REGISTER_CIRCUIT(rfedd);
  REGISTER_CIRCUIT(rlcg);
  REGISTER_CIRCUIT(spdeembed);
  REGISTER_CIRCUIT(spembed);
  REGISTER_CIRCUIT(spiralinductor);
  REGISTER_CIRCUIT(strafo);
  REGISTER_CIRCUIT(taperedline);
  REGISTER_CIRCUIT(thyristor);
  REGISTER_CIRCUIT(tline);
  REGISTER_CIRCUIT(tline4p);
  REGISTER_CIRCUIT(trafo);
  REGISTER_CIRCUIT(triac);
  REGISTER_CIRCUIT(tswitch);
  REGISTER_CIRCUIT(tunneldiode);
  REGISTER_CIRCUIT(twistedpair);
  REGISTER_CIRCUIT(vac);
  REGISTER_CIRCUIT(vam);
  REGISTER_CIRCUIT(vccs);
  REGISTER_CIRCUIT(vcvs);
  REGISTER_CIRCUIT(vdc);
  REGISTER_CIRCUIT(vexp);
  REGISTER_CIRCUIT(vfile);
  REGISTER_CIRCUIT(vnoise);
  REGISTER_CIRCUIT(vpm);
  REGISTER_CIRCUIT(vprobe);
  REGISTER_CIRCUIT(vpulse);
  REGISTER_CIRCUIT(vrect);
  REGISTER_CIRCUIT(vvnoise);
  REGISTER_CIRCUIT(wprobe);

  REGISTER_CIRCUIT(ecvs);

  // analyses
  REGISTER_ANALYSIS(acsolver);
  REGISTER_ANALYSIS(dcsolver);
  REGISTER_ANALYSIS(hbsolver);
  REGISTER_ANALYSIS(parasweep);
  REGISTER_ANALYSIS(spsolver);
  REGISTER_ANALYSIS(trsolver);
  REGISTER_ANALYSIS(trsolver);
}

} // namespace qucs
