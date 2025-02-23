/*
 * Copyright (C) 2003-2009 Stefan Jahn <stefan@lkcc.org>
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

#include <fstream>
#include <iostream>
#include <list>

#include "check_netlist.h"
#include "component.h"
#include "components.h"
#include "dataset.h"
#include "environment.h"
#include "exceptionstack.h"
#include "input.h"
#include "logging.h"
#include "module.h"
#include "net.h"

using namespace qucs;

int main(int argc, char **argv) {
  loginit();

  int ret = 0;

  char *infile = nullptr;
  char *outfile = nullptr;

  ::srand(::time(nullptr));

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      fprintf(stdout,
              "Usage: %s [OPTION]...\n\n"
              "  -h, --help     display this help and exit\n"
              "  -i FILENAME    use file as input netlist (default stdin)\n"
              "  -o FILENAME    use file as output dataset (default stdout)\n"
              "  -c, --check    check the input netlist and exit\n",
              argv[0]);
      return 0;
    } else if (!strcmp(argv[i], "-i")) {
      infile = argv[++i];
    } else if (!strcmp(argv[i], "-o")) {
      outfile = argv[++i];
      redirect_status_to_stdout();
    } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--check")) {
      netlist_check = 1;
    }
  }

  // create static modules
  module::registerModules();

  // create root environment
  const auto root = new environment(std::string("root"));

  // create netlist object and input
  const auto subnet = new net("subnet");
  const auto in = infile ? new input(infile) : new input();

  // pass root environment to netlist object and input
  subnet->setEnv(root);
  in->setEnv(root);

  // get input netlist
  if (in->netlist(subnet) != 0) {
    if (netlist_check) {
      logprint(LOG_STATUS, "checker notice, netlist check FAILED\n");
    }
    return -1;
  }
  if (netlist_check) {
    logprint(LOG_STATUS, "checker notice, netlist OK\n");
    return 0;
  }

  // attach a ground to the netlist
  const auto gnd = new ground();
  gnd->setNode(0, "gnd");
  gnd->setName("GND");
  subnet->insertCircuit(gnd);

  // analyse the netlist
  int err = 0;
  dataset *out = subnet->runAnalysis(err);
  ret |= err;

  // evaluate output dataset
  ret |= root->equationSolver(out);
  out->setFile(outfile);
  out->print();

  estack.print("uncaught");

  delete subnet;
  delete in;
  delete out;
  delete root;

  netlist_destroy_env();

  return ret;
}
