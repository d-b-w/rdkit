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
      if (subs[i].is_lone_pair) {
        // Lone pair has lowest priority (Z=0)
        priorities.emplace_back(0, i);
      } else {
        // Implicit H (Z=1, mass=1)
        priorities.emplace_back(1001, i);
      }
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
  // Implement CIP Rule 2: Compare paths incrementally
  // OPTIMIZATION: Build comparison keys incrementally by only processing NEW shell

  // Note: Substituents can have different numbers of shells (e.g., implicit H stops at shell 0)
  // Missing shells are treated as empty, which makes those subs lower priority

  // Append descriptors for the NEW shell to each substituent's comparison_key
  for (size_t i = 0; i < subs.size(); ++i) {
    // Skip if this substituent doesn't have this shell yet
    if (shell_depth >= subs[i].shells.size()) {
      continue;
    }

    const auto& shell = subs[i].shells[shell_depth];

    // Create atom descriptors with multiplicities for THIS shell only
    std::vector<int> shell_descriptors;
    shell_descriptors.reserve(shell.atoms.size());

    for (size_t j = 0; j < shell.atoms.size(); ++j) {
      const Atom* atom = shell.atoms[j];
      uint32_t mult = shell.multiplicities[j];

      int descriptor;
      if (atom == nullptr) {
        // Check if this is a lone pair or implicit H
        if (subs[i].is_lone_pair) {
          // Lone pair: Z=0 (lowest priority)
          descriptor = mult;
        } else {
          // Implicit H: Z=1, mass=1
          descriptor = 1001000 + mult;
        }
      } else {
        int z = atom->getAtomicNum();
        auto* pt = PeriodicTable::getTable();
        unsigned int mass = atom->getIsotope() ? atom->getIsotope() :
            static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z));

        // Encode as: Z * 1000000 + mass * 1000 + multiplicity
        // Multiplicity 0 means duplicate (ring closure)
        descriptor = z * 1000000 + static_cast<int>(mass) * 1000 + static_cast<int>(mult);
      }

      shell_descriptors.push_back(descriptor);
    }

    // Sort descriptors for this shell (descending = higher priority first)
    std::sort(shell_descriptors.rbegin(), shell_descriptors.rend());

    // Append to existing comparison key (INCREMENTAL!)
    subs[i].comparison_key.insert(subs[i].comparison_key.end(),
                                   shell_descriptors.begin(),
                                   shell_descriptors.end());
  }

  // Build sorting vector with references to existing keys
  std::vector<std::pair<const std::vector<int>*, size_t>> keys;
  keys.reserve(subs.size());

  for (size_t i = 0; i < subs.size(); ++i) {
    keys.emplace_back(&subs[i].comparison_key, i);
  }

  // Sort substituents by their keys (lexicographic, reverse for CIP priority)
  std::sort(keys.begin(), keys.end(),
            [](const auto& a, const auto& b) {
              return *a.first > *b.first;
            });

  // Check for ties
  for (size_t i = 1; i < keys.size(); ++i) {
    if (*keys[i].first == *keys[i-1].first) {
      return false;  // Still have ties, need more shells
    }
  }

  // All unique - assign ranks (rank 0 = highest priority)
  for (size_t rank = 0; rank < keys.size(); ++rank) {
    size_t sub_idx = keys[rank].second;
    subs[sub_idx].final_rank = static_cast<int>(rank);
  }

  return true;
}

namespace {

// Build comparison key for a substituent under a stereochemical hypothesis
// hypothesis: Descriptor::R or Descriptor::S - what to assume for unlabeled centers
std::vector<int> buildKeyWithHypothesis(const Substituent& sub,
                                        uint32_t max_shell,
                                        Descriptor hypothesis) {
  std::vector<int> key;

  for (uint32_t shell = 0; shell <= max_shell && shell < sub.shells.size(); ++shell) {
    const auto& shell_data = sub.shells[shell];

    std::vector<int> shell_descriptors;
    shell_descriptors.reserve(shell_data.atoms.size());

    for (size_t j = 0; j < shell_data.atoms.size(); ++j) {
      const Atom* atom = shell_data.atoms[j];
      uint32_t mult = shell_data.multiplicities[j];
      Descriptor stereo = shell_data.stereo_labels[j];

      int descriptor;
      if (atom == nullptr) {
        if (sub.is_lone_pair) {
          descriptor = mult;  // Z=0
        } else {
          descriptor = 1001000 + mult;  // Implicit H: Z=1
        }
      } else {
        int z = atom->getAtomicNum();
        auto* pt = PeriodicTable::getTable();
        unsigned int mass = atom->getIsotope() ? atom->getIsotope() :
            static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z));

        // Base descriptor: Z * 1000000 + mass * 1000 + multiplicity
        descriptor = z * 1000000 + static_cast<int>(mass) * 1000 + static_cast<int>(mult);

        // Add stereo contribution (Rule 5)
        // Encode stereo as additional high bits
        int stereo_value = 0;
        if (stereo != Descriptor::NONE) {
          // Use actual label if present
          if (stereo == Descriptor::R || stereo == Descriptor::r) {
            stereo_value = 2;  // R > S
          } else if (stereo == Descriptor::S || stereo == Descriptor::s) {
            stereo_value = 1;  // S < R
          }
          // E/Z also contributes
          else if (stereo == Descriptor::E) {
            stereo_value = 2;
          } else if (stereo == Descriptor::Z) {
            stereo_value = 1;
          }
        } else {
          // No label - use hypothesis for tetrahedral centers
          // We can detect tetrahedral by checking atom properties, but for simplicity,
          // apply hypothesis to all unlabeled centers
          if (hypothesis == Descriptor::R) {
            stereo_value = 2;
          } else if (hypothesis == Descriptor::S) {
            stereo_value = 1;
          }
        }

        // Add stereo as high-order bits (multiply by large number to keep it significant)
        descriptor += stereo_value * 1000000000;
      }

      shell_descriptors.push_back(descriptor);
    }

    // Sort shell descriptors (descending = higher priority first)
    std::sort(shell_descriptors.rbegin(), shell_descriptors.rend());

    // Append to key
    key.insert(key.end(), shell_descriptors.begin(), shell_descriptors.end());
  }

  return key;
}

}  // namespace

bool applyRule5(const ROMol& mol,
                std::vector<Substituent>& subs,
                uint32_t shell_depth,
                bool& is_pseudo) {
  if (subs.empty()) {
    return true;
  }

  // Build comparison keys under R hypothesis
  std::vector<std::vector<int>> keys_R;
  keys_R.reserve(subs.size());
  for (const auto& sub : subs) {
    keys_R.push_back(buildKeyWithHypothesis(sub, shell_depth, Descriptor::R));
  }

  // Build comparison keys under S hypothesis
  std::vector<std::vector<int>> keys_S;
  keys_S.reserve(subs.size());
  for (const auto& sub : subs) {
    keys_S.push_back(buildKeyWithHypothesis(sub, shell_depth, Descriptor::S));
  }

  // Build sorting indices with R hypothesis
  std::vector<std::pair<const std::vector<int>*, size_t>> sort_R;
  sort_R.reserve(subs.size());
  for (size_t i = 0; i < subs.size(); ++i) {
    sort_R.emplace_back(&keys_R[i], i);
  }

  std::sort(sort_R.begin(), sort_R.end(),
            [](const auto& a, const auto& b) {
              return *a.first > *b.first;
            });

  // Build sorting indices with S hypothesis
  std::vector<std::pair<const std::vector<int>*, size_t>> sort_S;
  sort_S.reserve(subs.size());
  for (size_t i = 0; i < subs.size(); ++i) {
    sort_S.emplace_back(&keys_S[i], i);
  }

  std::sort(sort_S.begin(), sort_S.end(),
            [](const auto& a, const auto& b) {
              return *a.first > *b.first;
            });

  // Check if R and S hypotheses give same ordering
  bool same_order = true;
  for (size_t i = 0; i < subs.size(); ++i) {
    if (sort_R[i].second != sort_S[i].second) {
      same_order = false;
      break;
    }
  }

  // If orderings differ, this is pseudo-asymmetric
  is_pseudo = !same_order;

  // Check for ties in R hypothesis ordering
  for (size_t i = 1; i < sort_R.size(); ++i) {
    if (*sort_R[i].first == *sort_R[i-1].first) {
      return false;  // Still have ties
    }
  }

  // All unique - assign ranks based on R hypothesis
  for (size_t rank = 0; rank < sort_R.size(); ++rank) {
    size_t sub_idx = sort_R[rank].second;
    subs[sub_idx].final_rank = static_cast<int>(rank);
  }

  return true;
}

bool checkPseudoAsymmetry(const std::vector<Substituent>& subs) {
  // This is now handled by applyRule5 setting the is_pseudo flag
  // This function is kept for compatibility but not used
  return false;
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
