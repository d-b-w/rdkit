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
#include <iostream>

namespace RDKit {
namespace NewCIPLabeler {

CenterRanking rankSubstituents(const ROMol& mol,
                               const Atom* center,
                               std::vector<Substituent>& subs,
                               uint32_t max_shells,
                               bool use_rule5) {
  if (subs.empty()) {
    return CenterRanking{true, false, {}};
  }

  // Fast path: try constitutional rules only (no shell expansion)
  if (tryConstitutionalRanking(mol, subs)) {
    return buildRankingResult(subs, false);
  }

  // Need deeper analysis - expand shells iteratively
  // Each substituent needs its own visited tracking
  std::vector<boost::dynamic_bitset<>> visited_sets;
  visited_sets.reserve(subs.size());
  for (size_t i = 0; i < subs.size(); ++i) {
    visited_sets.emplace_back(mol.getNumAtoms());
    visited_sets[i].set(center->getIdx());  // Don't backtrack to center
  }

  // Safety limits: prevent infinite loops on pathological molecules
  constexpr uint32_t HARD_MAX_SHELLS = 200;

  // If max_shells is 0, use reasonable default
  uint32_t max_iter = (max_shells == 0) ? HARD_MAX_SHELLS : std::min(max_shells, HARD_MAX_SHELLS);

  bool is_pseudo = false;

  for (uint32_t shell = 0; shell < max_iter; ++shell) {
    // Expand all substituents to this shell
    for (size_t i = 0; i < subs.size(); ++i) {
      if (subs[i].shells.size() <= shell) {
        subs[i].expandNextShell(mol, visited_sets[i]);
      }
    }

    // Apply CIP rules at this depth
    bool resolved = applyRankingRules(mol, subs, shell);
    if (resolved) {
      return buildRankingResult(subs, is_pseudo);
    }

    // If normal rules failed and Rule 5 is enabled, try stereochemical comparison
    if (use_rule5) {
      resolved = applyRule5(mol, subs, shell, is_pseudo);
      if (resolved) {
        return buildRankingResult(subs, is_pseudo);
      }
    }

    // Debug: warn if taking too long
    if (shell == 20 || shell == 50 || shell == 100 || shell == 150) {
      std::cerr << "WARNING: CIP ranking at center " << center->getIdx()
                << " reached shell " << shell << " without resolving\n";
      std::cerr << "  Center atom: " << center->getSymbol()
                << " with " << subs.size() << " substituents\n";
      for (size_t i = 0; i < subs.size(); ++i) {
        std::cerr << "    Sub[" << i << "]: ";
        if (subs[i].root_atom == nullptr) {
          std::cerr << (subs[i].is_lone_pair ? "LP" : "implH");
        } else {
          std::cerr << subs[i].root_atom->getSymbol() << " (idx=" << subs[i].root_atom->getIdx() << ")";
        }
        std::cerr << " - " << subs[i].shells.size() << " shells expanded\n";
      }
    }
  }

  // Exceeded max iterations - provide diagnostic info
  std::cerr << "ERROR: Max iterations exceeded at center " << center->getIdx()
            << " (" << center->getSymbol() << ")\n";
  std::cerr << "  This likely indicates a perfectly symmetric structure\n";
  std::cerr << "  or a bug in shell expansion. Substituents:\n";
  for (size_t i = 0; i < subs.size(); ++i) {
    std::cerr << "    [" << i << "]: ";
    if (subs[i].root_atom == nullptr) {
      std::cerr << (subs[i].is_lone_pair ? "LonePair" : "ImplicitH");
    } else {
      std::cerr << subs[i].root_atom->getSymbol() << "(idx=" << subs[i].root_atom->getIdx() << ")";
    }
    std::cerr << " - expanded to " << subs[i].shells.size() << " shells\n";
  }

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
