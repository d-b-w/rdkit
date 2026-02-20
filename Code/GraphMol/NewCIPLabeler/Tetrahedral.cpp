//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Tetrahedral.h"
#include "Ranker.h"
#include "Descriptor.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <algorithm>
#include <vector>

namespace RDKit {
namespace NewCIPLabeler {

namespace {

// Compute parity of permutation from reference order to CIP priority order
// Returns: 1 (odd permutation) or 2 (even permutation), 0 on error
int computeParity4(const std::vector<const Atom*>& cip_order,
                   const std::vector<const Atom*>& spatial_order) {
  if (cip_order.size() != 4 || spatial_order.size() != 4) {
    return 0;
  }

  // Build permutation array: for each spatial position, find its CIP rank
  std::vector<int> perm(4);
  for (size_t i = 0; i < 4; ++i) {
    // Find where spatial_order[i] appears in cip_order
    auto it = std::find(cip_order.begin(), cip_order.end(), spatial_order[i]);
    if (it == cip_order.end()) {
      return 0;  // Mismatch
    }
    perm[i] = static_cast<int>(std::distance(cip_order.begin(), it));
  }

  // Count inversions to determine parity
  int inversions = 0;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = i + 1; j < 4; ++j) {
      if (perm[i] > perm[j]) {
        inversions++;
      }
    }
  }

  // Odd inversions = odd permutation (1), even inversions = even permutation (2)
  return (inversions % 2 == 1) ? 1 : 2;
}

Descriptor computeTetrahedralDescriptor(const ROMol& mol,
                                        const Atom* center,
                                        const CenterRanking& ranking,
                                        const std::vector<Substituent>& subs) {
  if (!ranking.is_unique) {
    return Descriptor::NONE;
  }
  if (ranking.order.size() != 4) {
    return Descriptor::NONE;
  }

  // Get chiral tag
  auto tag = center->getChiralTag();
  if (tag != Atom::CHI_TETRAHEDRAL_CW && tag != Atom::CHI_TETRAHEDRAL_CCW) {
    return Descriptor::NONE;
  }

  // Build CIP priority order (high to low priority)
  // ranking.order contains substituent indices in decreasing priority (highest first)
  std::vector<const Atom*> cip_order;
  cip_order.reserve(4);
  for (const auto& idx : ranking.order) {
    cip_order.push_back(subs[idx].root_atom);
  }

  // Build spatial order (from bond iteration order)
  std::vector<const Atom*> spatial_order;
  spatial_order.reserve(4);
  for (const auto& bond : mol.atomBonds(center)) {
    spatial_order.push_back(bond->getOtherAtom(center));
  }
  // Add implicit H if needed
  if (spatial_order.size() == 3) {
    spatial_order.push_back(center);  // Use center as placeholder for implicit H
  }

  if (spatial_order.size() != 4) {
    return Descriptor::NONE;
  }

  // Compute parity
  int parity = computeParity4(cip_order, spatial_order);
  if (parity == 0) {
    return Descriptor::NONE;
  }

  // Adjust chiral tag based on parity
  auto config = tag;
  if (parity == 1) {  // Odd permutation - flip tag
    config = (tag == Atom::CHI_TETRAHEDRAL_CW) ?
             Atom::CHI_TETRAHEDRAL_CCW : Atom::CHI_TETRAHEDRAL_CW;
  }

  // Assign descriptor based on final configuration
  // CCW (@) → S, CW (@@) → R (matches old CIPLabeler)
  if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    return ranking.is_pseudo ? Descriptor::s : Descriptor::S;
  } else if (config == Atom::CHI_TETRAHEDRAL_CW) {
    return ranking.is_pseudo ? Descriptor::r : Descriptor::R;
  }

  return Descriptor::NONE;
}

}  // namespace

void labelTetrahedralCenter(ROMol& mol, Atom* center, uint32_t max_iters) {
  // Collect substituents
  std::vector<Substituent> subs;

  for (const auto& bond : mol.atomBonds(center)) {
    const Atom* nbr = bond->getOtherAtom(center);
    subs.emplace_back();
    subs.back().root_atom = nbr;
    subs.back().connecting_bond = bond;
  }

  // Add implicit H if present
  if (center->getTotalNumHs() > 0 && subs.size() == 3) {
    subs.emplace_back();
    subs.back().root_atom = nullptr;  // nullptr represents implicit H
    subs.back().connecting_bond = nullptr;
  }

  if (subs.size() != 4) {
    return;  // Not tetrahedral
  }

  // Rank substituents
  CenterRanking ranking = rankSubstituents(mol, center, subs, max_iters);

  if (!ranking.is_unique) {
    return;  // Cannot determine label
  }

  // Compute and assign descriptor
  Descriptor desc = computeTetrahedralDescriptor(mol, center, ranking, subs);
  if (desc != Descriptor::NONE) {
    std::string desc_str = to_string(desc);
    if (!desc_str.empty()) {
      center->setProp(common_properties::_CIPCode, desc_str);
    }
  }
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
