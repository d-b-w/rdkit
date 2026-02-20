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
  std::cerr << "[DEBUG]   computeTetrahedralDescriptor: ENTER" << std::endl;

  if (!ranking.is_unique) {
    std::cerr << "[DEBUG]     EXIT: ranking not unique" << std::endl;
    return Descriptor::NONE;
  }
  if (ranking.order.size() != 4) {
    std::cerr << "[DEBUG]     EXIT: order.size=" << ranking.order.size() << " (need 4)" << std::endl;
    return Descriptor::NONE;
  }

  // Get chiral tag
  auto tag = center->getChiralTag();
  std::cerr << "[DEBUG]     Chiral tag: " << static_cast<int>(tag) << std::endl;
  if (tag != Atom::CHI_TETRAHEDRAL_CW && tag != Atom::CHI_TETRAHEDRAL_CCW) {
    std::cerr << "[DEBUG]     EXIT: tag is not CW or CCW" << std::endl;
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

  std::cerr << "[DEBUG]     spatial_order.size=" << spatial_order.size() << std::endl;
  if (spatial_order.size() != 4) {
    std::cerr << "[DEBUG]     EXIT: spatial order not 4" << std::endl;
    return Descriptor::NONE;
  }

  // Compute parity
  int parity = computeParity4(cip_order, spatial_order);
  std::cerr << "[DEBUG]     Parity: " << parity << std::endl;
  if (parity == 0) {
    std::cerr << "[DEBUG]     EXIT: parity is 0" << std::endl;
    return Descriptor::NONE;
  }

  // Adjust chiral tag based on parity
  auto config = tag;
  if (parity == 1) {  // Odd permutation - flip tag
    config = (tag == Atom::CHI_TETRAHEDRAL_CW) ?
             Atom::CHI_TETRAHEDRAL_CCW : Atom::CHI_TETRAHEDRAL_CW;
  }

  std::cerr << "[DEBUG]     Original tag: " << static_cast<int>(tag)
            << ", config after parity: " << static_cast<int>(config) << std::endl;

  // Assign descriptor based on final configuration
  // CCW (@) → S, CW (@@) → R (matches old CIPLabeler)
  Descriptor result;
  if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    result = ranking.is_pseudo ? Descriptor::s : Descriptor::S;
  } else if (config == Atom::CHI_TETRAHEDRAL_CW) {
    result = ranking.is_pseudo ? Descriptor::r : Descriptor::R;
  } else {
    result = Descriptor::NONE;
  }

  std::cerr << "[DEBUG]     Returning descriptor: " << static_cast<int>(result) << std::endl;
  return result;
}

}  // namespace

void labelTetrahedralCenter(ROMol& mol, Atom* center, uint32_t max_iters) {
  std::cerr << "[DEBUG] labelTetrahedralCenter: ENTER for atom " << center->getIdx() << std::endl;
  std::cerr << "[DEBUG]   Chiral tag: " << static_cast<int>(center->getChiralTag()) << std::endl;

  // Check existing CIP code
  std::string existing_code;
  bool has_existing = center->getPropIfPresent(common_properties::_CIPCode, existing_code);
  std::cerr << "[DEBUG]   Existing _CIPCode: " << (has_existing ? existing_code : "NONE") << std::endl;

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

  std::cerr << "[DEBUG]   Collected " << subs.size() << " substituents" << std::endl;

  if (subs.size() != 4) {
    std::cerr << "[DEBUG]   EXIT: Not tetrahedral (need 4 substituents)" << std::endl;
    return;  // Not tetrahedral
  }

  // Rank substituents
  std::cerr << "[DEBUG]   Calling rankSubstituents..." << std::endl;
  CenterRanking ranking = rankSubstituents(mol, center, subs, max_iters);

  std::cerr << "[DEBUG]   Ranking result: is_unique=" << ranking.is_unique
            << ", order.size()=" << ranking.order.size() << std::endl;
  if (ranking.is_unique) {
    std::cerr << "[DEBUG]   Final ranks: ";
    for (size_t i = 0; i < subs.size(); ++i) {
      std::cerr << "sub[" << i << "].rank=" << subs[i].final_rank << " ";
    }
    std::cerr << std::endl;
  }

  if (!ranking.is_unique) {
    std::cerr << "[DEBUG]   EXIT: Ranking not unique" << std::endl;
    return;  // Cannot determine label
  }

  // Compute and assign descriptor
  std::cerr << "[DEBUG]   Calling computeTetrahedralDescriptor..." << std::endl;
  Descriptor desc = computeTetrahedralDescriptor(mol, center, ranking, subs);
  std::cerr << "[DEBUG]   Computed descriptor: " << static_cast<int>(desc)
            << " (" << to_string(desc) << ")" << std::endl;

  if (desc != Descriptor::NONE) {
    std::string desc_str = to_string(desc);
    std::cerr << "[DEBUG]   desc_str=\"" << desc_str << "\", empty=" << desc_str.empty() << std::endl;
    if (!desc_str.empty()) {
      std::cerr << "[DEBUG]   Setting _CIPCode to \"" << desc_str << "\"" << std::endl;
      center->setProp(common_properties::_CIPCode, desc_str);
      std::cerr << "[DEBUG]   setProp completed" << std::endl;
    } else {
      std::cerr << "[DEBUG]   EXIT: desc_str is empty!" << std::endl;
    }
  } else {
    std::cerr << "[DEBUG]   EXIT: Descriptor is NONE" << std::endl;
  }

  std::cerr << "[DEBUG] labelTetrahedralCenter: EXIT" << std::endl;
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
