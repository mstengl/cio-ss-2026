#pragma once
#include <objscip/objsepa.h>
#include <scip/type_retcode.h>

#include <optional>
#include <span>
#include <utility>
#include <vector>

class LiftedKnapsackSepa : public scip::ObjSepa {
 public:
  LiftedKnapsackSepa(SCIP* scip)
      : scip::ObjSepa(scip, "liftedknapsack", "Lifted Knapsack Separator",
                      99999, 1, 1.0, 0, false) {};
  SCIP_DECL_SEPAEXECLP(scip_execlp) override;

 private:
  void add_cut(SCIP* scip, SCIP_SEPA* sepa, std::span<int> indices,
               std::span<SCIP_Real> coefficients, SCIP_Real lhs, SCIP_Real rhs,
               bool forcecut, SCIP_Result* result);
  int cut_count = 0;
  std::optional<std::vector<SCIP_Row*>> knapsack_constraints = std::nullopt;

  /* Scratch buffers reused across rows and across calls; they only ever grow,
   * so the per-row work stays allocation-free after the first few rows. */
  std::vector<SCIP_Real> incumbent_solution;
  std::vector<int> indices;
  std::vector<int> items_taken;
  std::vector<int> items_not_taken;
  std::vector<std::pair<SCIP_Real, int>> cpuw_list;
  std::vector<SCIP_Real> coeffs;
  std::vector<int> original_indices;
};