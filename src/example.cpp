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

  auto env = new environment(std::string("root"));
  auto root = new net("subnet");
  root->setEnv(env);

  const auto in = new input("input.txt");
  in->setEnv(env);
  if (in->netlist(root) != 0) {
    if (netlist_check) {
      logprint(LOG_STATUS, "checker notice, netlist check FAILED\n");
    }
    return -1;
  }
  if (netlist_check) {
    logprint(LOG_STATUS, "checker notice, netlist OK\n");
    return 0;
  }

  auto *gnd = new ground();
  gnd->setNode(0, "gnd");
  gnd->setName("GND");
  root->insertCircuit(gnd);

  // analyse the netlist
  int err = 0;
  dataset *out = root->runAnalysis(err);
  ret |= err;

  // evaluate output dataset
  ret |= env->equationSolver(out);
  out->setFile("output.txt");
  out->print();

  estack.print("uncaught");

  delete in;
  delete out;
  delete root;
  delete env;

  netlist_destroy_env();

  return ret;
}
