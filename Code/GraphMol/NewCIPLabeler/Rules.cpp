//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Rules.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/PeriodicTable.h>
#include <algorithm>
#include <vector>
#include <map>

namespace RDKit {
namespace NewCIPLabeler {

namespace {

// Compare two atoms by CIP priority (higher atomic number/mass = higher priority)
// Returns: -1 if a < b, 0 if a == b, 1 if a > b
int compareAtomsPriority(const Atom* a, const Atom* b) {
  // Handle nullptr (implicit H)
  if (a == nullptr && b == nullptr) return 0;
  if (a == nullptr) return -1;  // H has lowest priority
  if (b == nullptr) return 1;

  // Rule 1a: Higher atomic number has priority
  int z_a = a->getAtomicNum();
  int z_b = b->getAtomicNum();
  if (z_a != z_b) {
    return (z_a > z_b) ? 1 : -1;
  }

  // Rule 1b: Higher isotope mass has priority
  auto* pt = PeriodicTable::getTable();
  unsigned int mass_a = a->getIsotope() ? a->getIsotope() :
      static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z_a));
  unsigned int mass_b = b->getIsotope() ? b->getIsotope() :
      static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z_b));

  if (mass_a != mass_b) {
    return (mass_a > mass_b) ? 1 : -1;
  }

  return 0;  // Equal priority
}

}  // namespace

bool tryConstitutionalRanking(const ROMol& mol,
                              std::vector<Substituent>& subs) {
  if (subs.empty()) return true;

  // Build sorting keys with original indices
  std::vector<std::pair<int, size_t>> priorities;
  priorities.reserve(subs.size());

  for (size_t i = 0; i < subs.size(); ++i) {
    // Use a simple priority score: atomic_num * 1000 + mass
    const Atom* atom = subs[i].root_atom;

    if (atom == nullptr) {
      // Implicit H
      priorities.emplace_back(1001, i);  // Z=1, mass=1
    } else {
      int z = atom->getAtomicNum();
      auto* pt = PeriodicTable::getTable();
      unsigned int mass = atom->getIsotope() ? atom->getIsotope() :
          static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z));
      int score = z * 1000 + static_cast<int>(mass);
      priorities.emplace_back(score, i);
    }
  }

  // Sort by priority (higher = higher priority, so reverse order)
  std::sort(priorities.begin(), priorities.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  // Check for ties
  for (size_t i = 1; i < priorities.size(); ++i) {
    if (priorities[i].first == priorities[i-1].first) {
      return false;  // Tie detected, need deeper analysis
    }
  }

  // All unique - assign ranks (rank 0 = highest priority)
  for (size_t rank = 0; rank < priorities.size(); ++rank) {
    size_t sub_idx = priorities[rank].second;
    subs[sub_idx].final_rank = static_cast<int>(rank);
  }

  return true;
}

bool applyRankingRules(const ROMol& mol,
                       std::vector<Substituent>& subs,
                       uint32_t shell_depth) {
  // For now, implement simple lexicographic comparison at the shell level
  // TODO: Full CIP rules implementation in Phase 2

  // If we don't have enough shells yet, can't resolve
  for (const auto& sub : subs) {
    if (sub.shells.size() <= shell_depth) {
      return false;
    }
  }

  // Build comparison keys for each substituent at this shell depth
  std::vector<std::pair<std::vector<int>, size_t>> keys;
  keys.reserve(subs.size());

  for (size_t i = 0; i < subs.size(); ++i) {
    const auto& shell = subs[i].shells[shell_depth];
    std::vector<int> key;

    // Create sorted list of (atomic_num, mass) pairs
    for (const auto* atom : shell.atoms) {
      if (atom == nullptr) {
        key.push_back(1001);  // H
      } else {
        int z = atom->getAtomicNum();
        auto* pt = PeriodicTable::getTable();
        unsigned int mass = atom->getIsotope() ? atom->getIsotope() :
            static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z));
        key.push_back(z * 1000 + static_cast<int>(mass));
      }
    }

    // Sort the key for lexicographic comparison
    std::sort(key.rbegin(), key.rend());  // Reverse sort (higher priority first)

    keys.emplace_back(std::move(key), i);
  }

  // Sort by keys
  std::sort(keys.begin(), keys.end(),
            [](const auto& a, const auto& b) {
              // Lexicographic comparison (reverse for priority)
              return a.first > b.first;
            });

  // Check for ties
  for (size_t i = 1; i < keys.size(); ++i) {
    if (keys[i].first == keys[i-1].first) {
      return false;  // Still have ties
    }
  }

  // All unique - assign ranks
  for (size_t rank = 0; rank < keys.size(); ++rank) {
    size_t sub_idx = keys[rank].second;
    subs[sub_idx].final_rank = static_cast<int>(rank);
  }

  return true;
}

bool checkPseudoAsymmetry(const std::vector<Substituent>& subs) {
  // Pseudo-asymmetry detection
  // This is a placeholder - full implementation in Phase 4
  return false;
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
