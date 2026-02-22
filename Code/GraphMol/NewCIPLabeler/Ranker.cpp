//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Ranker.h"
#include "Rules.h"
#include "NewCIPLabeler.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>

namespace RDKit {
namespace NewCIPLabeler {

CenterRanking rankSubstituents(const ROMol& mol,
                               const Atom* center,
                               std::vector<Substituent>& subs,
                               uint32_t max_shells) {
  if (subs.empty()) {
    return CenterRanking{true, false, {}};
  }

  // Fast path: try constitutional rules only (no shell expansion)
  if (tryConstitutionalRanking(mol, subs)) {
    return buildRankingResult(subs, false);
  }

  // Need deeper analysis - expand shells iteratively
  boost::dynamic_bitset<> visited(mol.getNumAtoms());
  visited.set(center->getIdx());  // Don't backtrack to center

  // Safety limits: prevent infinite loops on pathological molecules
  constexpr uint32_t HARD_MAX_SHELLS = 100;

  // If max_shells is 0, use reasonable default
  uint32_t max_iter = (max_shells == 0) ? HARD_MAX_SHELLS : std::min(max_shells, HARD_MAX_SHELLS);

  for (uint32_t shell = 0; shell < max_iter; ++shell) {
    // Expand all substituents to this shell
    for (auto& sub : subs) {
      if (sub.shells.size() <= shell) {
        sub.expandNextShell(mol, visited);
      }
    }

    // Apply CIP rules at this depth
    bool resolved = applyRankingRules(mol, subs, shell);
    if (resolved) {
      bool is_pseudo = checkPseudoAsymmetry(subs);
      return buildRankingResult(subs, is_pseudo);
    }
  }

  // Exceeded max iterations
  throw MaxIterationsExceeded();
}

CenterRanking buildRankingResult(const std::vector<Substituent>& subs,
                                 bool is_pseudo) {
  CenterRanking result;
  result.is_pseudo = is_pseudo;

  // Check if all substituents are uniquely ranked
  std::map<int, int> rank_counts;
  for (const auto& sub : subs) {
    if (sub.final_rank < 0) {
      result.is_unique = false;
      return result;
    }
    rank_counts[sub.final_rank]++;
  }

  // Check for duplicate ranks
  for (const auto& [rank, count] : rank_counts) {
    if (count > 1) {
      result.is_unique = false;
      return result;
    }
  }

  result.is_unique = true;

  // Build order vector (indices sorted by rank)
  // Note: rank 0 = highest priority, so sorting ascending gives highest priority first
  std::vector<std::pair<int, size_t>> rank_idx;
  for (size_t i = 0; i < subs.size(); ++i) {
    rank_idx.emplace_back(subs[i].final_rank, i);
  }

  // Sort by rank (ascending = rank 0 first = highest priority first)
  std::sort(rank_idx.begin(), rank_idx.end());

  result.order.reserve(rank_idx.size());
  for (const auto& [rank, idx] : rank_idx) {
    result.order.push_back(idx);
  }

  return result;
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
