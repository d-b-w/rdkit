//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <vector>
#include <cstdint>
#include <RDGeneral/export.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "Descriptor.h"

namespace RDKit {

class ROMol;
class Atom;
class Bond;

namespace NewCIPLabeler {

// Single shell of atoms at distance N from root
struct AtomShell {
  uint32_t distance;
  std::vector<const Atom*> atoms;
  std::vector<uint32_t> multiplicities;  // Bond order duplicates
  std::vector<Descriptor> stereo_labels;  // R/S/E/Z labels for Rule 5
};

// One substituent extending from a stereocenter
struct Substituent {
  const Atom* root_atom;           // Immediate neighbor (nullptr for implicit H/lone pair)
  const Bond* connecting_bond;     // Bond to center (nullptr for implicit H/lone pair)
  std::vector<AtomShell> shells;   // Expanded shells (lazy)
  std::vector<int> comparison_key; // Incrementally built CIP comparison key
  int final_rank = -1;             // Assigned rank (-1 if not yet ranked)
  bool is_lone_pair = false;       // True if this is a lone pair (lower priority than H)

  // Expand one more shell
  void expandNextShell(const ROMol& mol, boost::dynamic_bitset<>& visited);
};

// Ranking result for a center
struct CenterRanking {
  bool is_unique;             // All substituents uniquely ranked?
  bool is_pseudo;             // Use r/s instead of R/S?
  std::vector<size_t> order;  // Sorted indices (low to high priority)
};

}  // namespace NewCIPLabeler
}  // namespace RDKit
