#include "heur_zirounding.hpp"

#include <scip/pub_lp.h>
#include <scip/pub_var.h>
#include <scip/scip_branch.h>
#include <scip/scip_lp.h>
#include <scip/scip_numerics.h>
#include <scip/scip_solvingstats.h>
#include <scip/type_heur.h>
#include <scip/type_result.h>
#include <scip/type_retcode.h>
#include <scip/type_var.h>

#include <cmath>
#include <cstddef>
#include <latch>
#include <memory>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "utils.hpp"

static constexpr SCIP_Real min_shift =
    1e-4;  // Minimum change to the incumbent for a move to be considered
static constexpr SCIP_Real min_slack =
    1e-5;  // Minimum slack that a constraint should have before move can be
           // generated from it
static constexpr SCIP_Real max_iter = 1'000'000;

/* Do not change start */
struct SolDeleter {
  SCIP* scip;
  SolDeleter(SCIP* scip) : scip(scip) {};
  void operator()(SCIP_Sol* p) { SCIPfreeSol(scip, &p); };
};
using SolPtr = std::unique_ptr<SCIP_SOL, SolDeleter>;
bool add_solution_helper(SCIP* scip, SCIP_HEUR* heur,
                         std::span<SCIP_COL* const> lp_cols,
                         std::span<SCIP_Real> incumbent) {
  SolDeleter deleter(scip);
  SolPtr sol(nullptr, deleter);
  CALL_CHECK(SCIPcreateSol(scip, std::out_ptr(sol), heur));
  for (auto [col, val] : std::ranges::zip_view(lp_cols, incumbent)) {
    CALL_CHECK(SCIPsetSolVal(scip, sol.get(), SCIPcolGetVar(col), val));
  }
  unsigned int stored;
  CALL_CHECK(
      SCIPtrySol(scip, sol.get(), FALSE, TRUE, TRUE, TRUE, TRUE, &stored));
  return stored;
}
/* Do not change end*/
bool is_any_fractional_remaining(SCIP* scip, std::span<SCIP_COL* const> lp_cols,
                                 std::span<SCIP_Real> incumbent) {
  // HINT use SCIPvarIsIntegral ad SCIPcolGetVar to check if any incumbent is
  // still fractional, use SCIPisFeasIntegral to check if value is within
  // integer tolerances
  for (auto [col, val] : std::ranges::zip_view(lp_cols, incumbent)) {
    if (!SCIPvarIsIntegral(SCIPcolGetVar(col))) continue;
    if (!SCIPisFeasIntegral(scip, val)) return true;
  }
  return false;
}

struct FractionalVarIdx {
  size_t index;
};
std::vector<FractionalVarIdx> collect_fractional_remaining(
    SCIP* scip, std::span<SCIP_COL* const> lp_cols,
    std::span<SCIP_Real> incumbent) {
  // HINT use SCIPvarIsIntegral ad SCIPcolGetVar to check if any incumbent is
  // still fractional, use SCIPisFeasIntegral to check if value is within
  // integer tolerances
  std::vector<FractionalVarIdx> fractional_vals;
  for (auto [col, val, i] :
       std::ranges::zip_view(lp_cols, incumbent,
                             std::ranges::views::iota(0U, incumbent.size()))) {
    if (!SCIPvarIsIntegral(SCIPcolGetVar(col))) continue;
    if (!SCIPisFeasIntegral(scip, val))
      fractional_vals.push_back(FractionalVarIdx{i});
  }
  return fractional_vals;
}

SCIP_Real compute_incumbent_activity(SCIP_Row* row,
                                     std::span<SCIP_Real> incumbent) {
  // Compute incumbent activity for the current row
  // HINT create a coeff span using SCIProwGetVals, SCIProwGetNNonz and a cols
  // span using SCIProwGetCols, and SCIProwGetNNonz. Finally use SCIPcolGetLPPos
  // row activity = row constant + sum over nonzero coeff*value. The constant
  // term can be obtianed via SCIProwGetConstant
  auto nnonz = static_cast<size_t>(SCIProwGetNNonz(row));
  std::span<SCIP_Real const> coeffs{SCIProwGetVals(row), nnonz};
  std::span<SCIP_COL* const> row_cols{SCIProwGetCols(row), nnonz};

  SCIP_Real activity = SCIProwGetConstant(row);
  for (auto [coeff, col] : std::ranges::zip_view(coeffs, row_cols)) {
    int lppos = SCIPcolGetLPPos(col);
    if (lppos < 0) continue;  // column not in the current LP
    activity += coeff * incumbent[static_cast<size_t>(lppos)];
  }
  return activity;
}

SCIP_DECL_HEUREXEC(ZIRoundHeur::scip_exec) {
  /* Dont touch start from here */
  *result = SCIP_DIDNOTRUN;
  assert(SCIPhasCurrentNodeLP(scip));
  assert(SCIPgetNLPs(scip) != nlpsolve);

  // Do not call heuristic if node was already detected to be infeasible
  if (nodeinfeasible) return SCIP_OKAY;

  // Only call heuristic if an optimal LP solution is at hand
  if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) return SCIP_OKAY;

  // Only call heuristic if LP objective value is smaller than the cutoff bound
  if (SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
    return SCIP_OKAY;

  // If LP Solution satisfy integrality constraint then skip the heuristic
  if (SCIPgetNLPBranchCands(scip) == 0) return SCIP_OKAY;

  // Make sure we are at a new LP solution before heuristic is called
  nlpsolve = SCIPgetNLPs(scip);

  *result = SCIP_DIDNOTFIND;
  /* Don't touch end here */

  /* Start editing from here onwards*/
  // Get LP informations, construct spans from SCIPgetLPRows, SCIPgetNLPRows,
  // SCIPgetLPCols, SCIPgetNLPCols
  std::span<SCIP_COL*> lp_cols{SCIPgetLPCols(scip),
                               static_cast<size_t>(SCIPgetNLPCols(scip))};
  std::span<SCIP_ROW*> lp_rows{SCIPgetLPRows(scip),
                               static_cast<size_t>(SCIPgetNLPRows(scip))};

  // Read the incumbent LP solution from SCIP see SCIPcolGetPrimsol
  std::vector<SCIP_Real> incumbent(lp_cols.size());
  std::ranges::transform(lp_cols, incumbent.begin(),
                         [](SCIP_COL* c) { return SCIPcolGetPrimsol(c); });
  // Loop while there is any fractional variable remaining in the incumbent
  // solution and iter have not reached max iter
  int iter = 0;
  auto fractional_vals = collect_fractional_remaining(scip, lp_cols, incumbent);
  bool shift_found = true;
  // Cache of the row activities under the current incumbent, indexed by LP row
  // position. Computed once here and kept in sync incrementally whenever a
  // shift is committed to the incumbent.
  std::vector<SCIP_Real> activities(lp_rows.size());
  std::ranges::transform(lp_rows, activities.begin(), [&](SCIP_ROW* r) {
    return compute_incumbent_activity(r, incumbent);
  });
  std::vector<size_t> to_remove_frac_indices;
  while (iter < max_iter && !fractional_vals.empty() && shift_found) {
    shift_found = false;
    for (const auto _i : std::ranges::views::iota(0U, fractional_vals.size())) {
      const auto var_idx = fractional_vals[_i];
      SCIP_COL* col = lp_cols[var_idx.index];
      SCIP_Real x = incumbent[var_idx.index];

      // start with the column's own bounds
      SCIP_Real var_ub = SCIPcolGetUb(col);
      SCIP_Real var_lb = SCIPcolGetLb(col);
      // rows in which this variable has a nonzero coefficient
      auto col_nnonz = static_cast<size_t>(SCIPcolGetNNonz(col));
      std::span<SCIP_Real const> var_coeffs{SCIPcolGetVals(col), col_nnonz};
      std::span<SCIP_ROW* const> rows_var_appears{SCIPcolGetRows(col),
                                                  col_nnonz};
      for (auto [a, row] :
           std::ranges::zip_view(var_coeffs, rows_var_appears)) {
        if (!SCIProwIsInLP(row)) continue;
        if (SCIPisZero(scip, a)) continue;

        SCIP_Real activity =
            activities[static_cast<size_t>(SCIProwGetLPPos(row))];
        SCIP_Real lhs = SCIProwGetLhs(row);
        SCIP_Real rhs = SCIProwGetRhs(row);

        // slack to each side; infinite side => no restriction
        SCIP_Real up_slack =
            SCIPisInfinity(scip, rhs) ? SCIPinfinity(scip) : rhs - activity;
        SCIP_Real down_slack =
            SCIPisInfinity(scip, -lhs) ? SCIPinfinity(scip) : activity - lhs;

        if (a > 0.0) {
          // increasing x raises activity -> limited by rhs; decreasing -> by
          // lhs
          if (!SCIPisInfinity(scip, up_slack))
            var_ub = MIN(var_ub, x + up_slack / a);
          if (!SCIPisInfinity(scip, down_slack))
            var_lb = MAX(var_lb, x - down_slack / a);
        } else {
          // negative coefficient: directions swap
          if (!SCIPisInfinity(scip, down_slack))
            var_ub = MIN(var_ub, x - down_slack / a);
          if (!SCIPisInfinity(scip, up_slack))
            var_lb = MAX(var_lb, x + up_slack / a);
        }
        if (SCIPisLE(scip, var_ub - var_lb, 0.0)) break;  // frozen, stop early
      }

      // check if the next integer up or down is reachable without breaking
      // anything use SCIPfeasCeil/Floor to get next integer up or down here,
      // since we don't need to ceil if the value is solution is already near
      // enough to integral To do LE and GE with tolerances use SCIPisLE and
      // SCIPisGE To do LT and GT with tolerances use SCIPisLT and SCIPisGT
      const SCIP_Real ceil_x = SCIPfeasCeil(scip, x);
      const SCIP_Real floor_x = SCIPfeasFloor(scip, x);
      const bool up_roundable = SCIPisLE(scip, ceil_x, var_ub);
      const bool down_roundable = SCIPisGE(scip, floor_x, var_lb);

      // SCIP minimises internally, so a positive objective coefficient favours
      // decreasing the variable and a negative one favours increasing it
      const SCIP_Real obj = SCIPcolGetObj(col);

      // fractionality of a value: distance to the nearest integer
      auto zi = [&](SCIP_Real v) {
        return MIN(v - SCIPfeasFloor(scip, v), SCIPfeasCeil(scip, v) - v);
      };

      // If var is up and down roundable pick objective improving direction
      // If var is only up or only down roundable pick roundable direction
      // If var is neither roundable determine whether staying, going to the
      // upper bound, going to the lower bound give the best fractionality
      // reduction Tie can be broken by selecting the objective improving
      // direction
      SCIP_Real new_val = x;
      if (up_roundable && down_roundable) {
        if (SCIPisPositive(scip, obj)) {
          new_val = floor_x;
        } else if (SCIPisNegative(scip, obj)) {
          new_val = ceil_x;
        } else {
          new_val = SCIPisLE(scip, ceil_x - x, x - floor_x) ? ceil_x : floor_x;
        }
      } else if (up_roundable) {
        new_val = ceil_x;
      } else if (down_roundable) {
        new_val = floor_x;
      } else {
        // Neither integer is reachable: move to whichever bound reduces the
        // fractionality most, staying put if neither of them helps
        SCIP_Real best_zi = zi(x);
        auto consider = [&](SCIP_Real cand) {
          if (SCIPisInfinity(scip, REALABS(cand))) return;
          const SCIP_Real cand_zi = zi(cand);
          const SCIP_Real delta_obj = obj * (cand - new_val);
          if (SCIPisLT(scip, cand_zi, best_zi) ||
              (SCIPisEQ(scip, cand_zi, best_zi) &&
               SCIPisNegative(scip, delta_obj))) {
            new_val = cand;
            best_zi = cand_zi;
          }
        };
        consider(var_ub);
        consider(var_lb);
      }

      // Check if we find any shift
      const SCIP_Real delta = new_val - x;
      if (REALABS(delta) >= min_shift) {
        incumbent[var_idx.index] = new_val;
        shift_found = true;

        // keep the cached row activities in sync with the incumbent
        for (auto [a, row] :
             std::ranges::zip_view(var_coeffs, rows_var_appears)) {
          if (!SCIProwIsInLP(row)) continue;
          activities[static_cast<size_t>(SCIProwGetLPPos(row))] += a * delta;
        }

        // drop the variable from the candidate list once it became integral
        if (SCIPisFeasIntegral(scip, new_val)) {
          to_remove_frac_indices.push_back(_i);
        }
      }

      // Increment iter_count and check if it is less than max_iter
      iter++;
      if (iter >= max_iter) break;
    }
    while (!to_remove_frac_indices.empty()) {
      const auto last_i = to_remove_frac_indices.back();
      to_remove_frac_indices.pop_back();
      std::swap(fractional_vals.at(last_i), fractional_vals.back());
      fractional_vals.pop_back();
    }
    // check if var is fractional see is_fractional_remaining hint
    // compute how far up / down this variable can be rounded
    // start with var_ub and var_lb via SCIPcolGetUb/ SCIPcolGetLb

    // for each row containing this variable, recompute the row's current
    // activity and see how it limits this variable's value
    // Hint use SCIPcolGetVals and SCIPcolGetNNonz to get the coefficient of
    // the columns on rows that are non zero and use SCIPcolGetRows to get
    // the rows for which the variable coefficients is nonzero

    // Loop over the rows update var_ub and var_lb
    // Hint 1: during the for loop over the rows check if row is in LP using
    // if (!SCIProwIsInLP(row)) continue;
    // Hint 2: Use SCIProwGetLhs and SCIProwGetRhs to get left and right and
    // side remember to check if they are + or - infity using
    // SCIPisInfinity, note it applies if l is - SCIPinfinity(scip.get())
    // then SCIPisInifinity(-l) is true
    // Hint 3: SCIP provide with MIN and MAX function
    // Determine update rule for var_ub and var_lb
  }

  // Add Sol if any is found
  // Use
  // bool stored = add_solution_helper(scip, heur, lp_cols, incumbent);
  // if(stored) *result = SCIP_FOUNDSOL;
  if (fractional_vals.empty()) {
    bool stored = add_solution_helper(scip, heur, lp_cols, incumbent);
    if (stored) {
      *result = SCIP_FOUNDSOL;
    }
  }

  return SCIP_OKAY;
}
