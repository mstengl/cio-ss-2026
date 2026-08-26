#include <objscip/objscip.h>
#include <objscip/objscipdefplugins.h>
#include <objscip/objsepa.h>
#include <scip/dialog_default.h>
#include <scip/scip_solvingstats.h>

#include <boost/program_options.hpp>
#include <cassert>
#include <format>
#include <iostream>
#include <memory>
#include <string>

#include "heur_zirounding.hpp"
#include "sepa_liftedknapsack.hpp"
#include "utils.hpp"

void scip_include_minimal_set(SCIP* scip) {
  CALL_CHECK(SCIPincludeDialogDefaultBasic(scip));
  CALL_CHECK(SCIPincludeConshdlrLinear(scip));
  CALL_CHECK(SCIPincludeConshdlrIntegral(scip));
  CALL_CHECK(SCIPincludeReaderMps(scip));
  CALL_CHECK(SCIPincludeNodeselBfs(scip));
  CALL_CHECK(SCIPincludeBranchruleRandom(scip));
  CALL_CHECK(SCIPincludeEventHdlrSolvingphase(scip));
  CALL_CHECK(SCIPincludeComprLargestrepr(scip));
  CALL_CHECK(SCIPincludeComprWeakcompr(scip));
  CALL_CHECK(SCIPincludeCutselHybrid(scip));
  CALL_CHECK(SCIPincludeDispDefault(scip));
  CALL_CHECK(SCIPincludeTableDefault(scip));
}

struct SCIPDeleter {
  void operator()(SCIP* scip) const { SCIPfree(&scip); }
};
using SCIPPtr = std::unique_ptr<SCIP, SCIPDeleter>;

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "Display this help menu")(
      "input_path,i", po::value<std::string>()->required(),
      "Path to the input MPS file")(
      "settings_path,s",
      po::value<std::string>()->default_value("./default.set")->required(),
      "A path to a .set file to be used by scip. Default: default.set");

  po::positional_options_description pos;
  pos.add("input_path", 1);

  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
      vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  SCIPPtr scip;
  CALL_CHECK(SCIPcreate(std::out_ptr(scip)));
  scip_include_minimal_set(scip.get());

  CALL_CHECK(SCIPincludeObjHeur(scip.get(), new ZIRoundHeur(scip.get()), TRUE))
  CALL_CHECK(
      SCIPincludeObjSepa(scip.get(), new LiftedKnapsackSepa(scip.get()), TRUE))
  CALL_CHECK(SCIPreadProb(scip.get(),
                          vm["input_path"].as<std::string>().c_str(), nullptr));

  CALL_CHECK(SCIPreadParams(scip.get(),
                            vm["settings_path"].as<std::string>().c_str()));

  CALL_CHECK(SCIPsolve(scip.get()));

  CALL_CHECK(SCIPprintStatistics(scip.get(), nullptr));

  // Pretty printing of some stuff we need for benchmark
  auto solved = (SCIPgetStatus(scip.get()) == SCIP_STATUS_OPTIMAL);
  std::cout << std::format("Solution Status: {}\n",
                           (solved) ? "solved" : "not solved");
  std::cout << std::format("Time: {}\n", SCIPgetSolvingTime(scip.get()));
  std::cout << std::format("Nodes: {}\n", SCIPgetNTotalNodes(scip.get()));
  std::cout << std::format("PrimalDualIntegral: {}\n",
                           SCIPgetPrimalDualIntegral(scip.get()));
}
