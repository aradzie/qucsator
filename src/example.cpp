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
  int ret = 0;

  loginit();

  // create static modules
  module::registerModules();

  // create root environment
  environment *root = new environment(std::string("root"));

  // create netlist object and input
  net *subnet = new net("subnet");
  input *in = new input("input.txt");

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
  out->setFile("output.txt");
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
