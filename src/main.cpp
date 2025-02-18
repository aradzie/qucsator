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

#include "config.h"

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

/*! \todo replace environment name root by "/" in order to be filesystem compatible */
int main(int argc, char **argv) {
  loginit();

  char *infile = nullptr;
  char *outfile = nullptr;
  char *projPath = nullptr;
  int listing = 0;
  int ret = 0;
  int dynamicLoad = 0;
  std::list<std::string> vamodules;

  ::srand(::time(nullptr));

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      fprintf(stdout,
              "Usage: %s [OPTION]...\n\n"
              "  -h, --help     display this help and exit\n"
              "  -i FILENAME    use file as input netlist (default stdin)\n"
              "  -o FILENAME    use file as output dataset (default stdout)\n"
              "  -c, --check    check the input netlist and exit\n"
#if DEBUG
              "  -l, --listing  emit C-code for available definitions\n"
#endif
              "  -p, --path     project path (or location of dynamic modules)\n"
              "  -m, --module   list of dynamic loaded modules (base names separated by space)\n",
              argv[0]);
      return 0;
    } else if (!strcmp(argv[i], "-i")) {
      infile = argv[++i];
    } else if (!strcmp(argv[i], "-o")) {
      outfile = argv[++i];
      redirect_status_to_stdout();
    } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--check")) {
      netlist_check = 1;
    } else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--listing")) {
      listing = 1;
    }
    // \todo remove command arguments
    else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--path")) {
      projPath = argv[++i];
    } else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--module")) {
      dynamicLoad = 1;
    } else {
      if (dynamicLoad) {
        vamodules.push_back(argv[i]);
      }
    }
  }

  // create static modules
  module::registerModules();

#ifndef NOMODULEPRINT
  // emit C-code for available definitions if requested and exit
  if (listing) {
    module::print();
    return ret;
  }
#endif

  // look for dynamic libs, load and register them
  // \todo, keep this way of loading or keep only annotated netlist?
  if (dynamicLoad) {
    module::registerDynamicModules(projPath, vamodules);
  } else {
    // no argument, look into netlist
    std::string sLine = "";
    std::ifstream file;

    std::string projPathNet = "";

    file.open(infile);

    while (!file.eof()) {
      getline(file, sLine);

      if (sLine.find("--path") != std::string::npos) {
        std::cout << sLine << std::endl;

        size_t pos = 0;
        pos = sLine.find("=");
        sLine.erase(0, pos + 1);
        std::cout << sLine << std::endl;

        // projPath =  const_cast<char*>(sLine.c_str());
        // projPath = (char*)sLine.c_str();
        // std::cout << "inside" << projPath << std::endl;

        projPathNet = sLine;
      }

      if (sLine.find("--module") != std::string::npos) {
        std::cout << sLine << std::endl;

        size_t pos = 0;
        pos = sLine.find("=");
        sLine.erase(0, pos + 1);
        // std::cout << sLine << std::endl;

        // put module names into list
        std::istringstream ss(sLine);
        std::string token;

        while (std::getline(ss, token, ' ')) {
          std::cout << token << '\n';

          vamodules.push_back(token);
        }
      }
    }

    module::registerDynamicModules((char *)projPathNet.c_str(), vamodules);
    file.close();
  }

  // create root environment
  environment *root = new environment(std::string("root"));

  // create netlist object and input
  net *subnet = new net("subnet");
  input *in = infile ? new input(infile) : new input();

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
  circuit *gnd = new ground();
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

  // delete static modules and dynamic modules
  module::unregisterModules();

  // close all the dynamic libs if any opened
  module::closeDynamicLibs();

  netlist_destroy_env();
  return ret;
}
