//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DoubleBond.h"
#include "Ranker.h"
#include "Descriptor.h"
#include "Substituent.h"
#include "NewCIPLabeler.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <RDGeneral/types.h>
#include <algorithm>
#include <vector>

namespace RDKit {
namespace NewCIPLabeler {

namespace {

// Determine E/Z based on CIP rankings and bond stereo
Descriptor computeDoubleBondDescriptor(const ROMol& mol,
                                       const Bond* bond,
                                       const std::vector<unsigned int>& ranked_atoms) {
  // ranked_atoms[0] = highest priority on begin atom side
  // ranked_atoms[1] = highest priority on end atom side

  if (ranked_atoms.size() != 2) {
    return Descriptor::NONE;
  }

  // Get bond stereo (STEREOCIS/STEREOTRANS or STEREOE/STEREOZ)
  auto stereo = bond->getStereo();

  // Map bond stereo to E/Z based on which atoms have highest priority
  // If stereo atoms match ranked atoms order → use stereo directly
  // Otherwise need to flip

  auto stereo_atoms = bond->getStereoAtoms();
  if (stereo_atoms.size() != 2) {
    return Descriptor::NONE;
  }

  bool same_order = (stereo_atoms[0] == static_cast<int>(ranked_atoms[0]) &&
                     stereo_atoms[1] == static_cast<int>(ranked_atoms[1]));

  // STEREOTRANS/STEREOE means reference atoms are trans
  // STEREOCIS/STEREOZ means reference atoms are cis
  if (stereo == Bond::STEREOTRANS || stereo == Bond::STEREOE) {
    // Reference atoms are trans
    // If ranked atoms are same as reference → E (trans)
    // If ranked atoms are opposite → Z (cis)
    return same_order ? Descriptor::E : Descriptor::Z;
  } else if (stereo == Bond::STEREOCIS || stereo == Bond::STEREOZ) {
    // Reference atoms are cis
    // If ranked atoms are same as reference → Z (cis)
    // If ranked atoms are opposite → E (trans)
    return same_order ? Descriptor::Z : Descriptor::E;
  }

  return Descriptor::NONE;
}

}  // namespace

void labelDoubleBond(ROMol& mol, Bond* bond, uint32_t max_iters) {
  // Must be a double bond
  if (bond->getBondType() != Bond::DOUBLE) {
    return;
  }

  // Must have stereo information
  auto stereo = bond->getStereo();
  if (stereo != Bond::STEREOCIS && stereo != Bond::STEREOTRANS &&
      stereo != Bond::STEREOE && stereo != Bond::STEREOZ) {
    return;
  }

  Atom* begin_atom = bond->getBeginAtom();
  Atom* end_atom = bond->getEndAtom();

  // Collect substituents for begin atom (excluding end atom)
  std::vector<Substituent> begin_subs;
  for (const auto& nbr_bond : mol.atomBonds(begin_atom)) {
    const Atom* nbr = nbr_bond->getOtherAtom(begin_atom);
    if (nbr == end_atom) continue;  // Skip the double bond partner

    begin_subs.emplace_back();
    begin_subs.back().root_atom = nbr;
    begin_subs.back().connecting_bond = nbr_bond;
  }

  // Add implicit H if present
  if (begin_atom->getTotalNumHs() > 0 && begin_subs.size() == 1) {
    begin_subs.emplace_back();
    begin_subs.back().root_atom = nullptr;  // Implicit H
    begin_subs.back().connecting_bond = nullptr;
  }

  // Collect substituents for end atom (excluding begin atom)
  std::vector<Substituent> end_subs;
  for (const auto& nbr_bond : mol.atomBonds(end_atom)) {
    const Atom* nbr = nbr_bond->getOtherAtom(end_atom);
    if (nbr == begin_atom) continue;  // Skip the double bond partner

    end_subs.emplace_back();
    end_subs.back().root_atom = nbr;
    end_subs.back().connecting_bond = nbr_bond;
  }

  // Add implicit H if present
  if (end_atom->getTotalNumHs() > 0 && end_subs.size() == 1) {
    end_subs.emplace_back();
    end_subs.back().root_atom = nullptr;  // Implicit H
    end_subs.back().connecting_bond = nullptr;
  }

  // Need exactly 2 substituents on each side for E/Z
  if (begin_subs.size() != 2 || end_subs.size() != 2) {
    return;
  }

  // Rank substituents on begin atom
  CenterRanking begin_ranking;
  try {
    begin_ranking = rankSubstituents(mol, begin_atom, begin_subs, max_iters);
  } catch (const MaxIterationsExceeded& e) {
    // Hit max iterations - likely a perfectly symmetric structure
    return;
  }
  if (!begin_ranking.is_unique) {
    return;  // Cannot determine label
  }

  // Rank substituents on end atom
  CenterRanking end_ranking;
  try {
    end_ranking = rankSubstituents(mol, end_atom, end_subs, max_iters);
  } catch (const MaxIterationsExceeded& e) {
    // Hit max iterations - likely a perfectly symmetric structure
    return;
  }
  if (!end_ranking.is_unique) {
    return;  // Cannot determine label
  }

  // Get highest priority atom indices
  // ranking.order[0] is index of highest priority substituent
  const Atom* begin_highest = begin_subs[begin_ranking.order[0]].root_atom;
  const Atom* end_highest = end_subs[end_ranking.order[0]].root_atom;

  // nullptr represents implicit H - need to get actual stereo atom from RDKit
  std::vector<unsigned int> ranked_atoms(2);

  // For begin atom side
  if (begin_highest == nullptr) {
    // Implicit H - RDKit uses the double bond partner as placeholder
    ranked_atoms[0] = end_atom->getIdx();
  } else {
    ranked_atoms[0] = begin_highest->getIdx();
  }

  // For end atom side
  if (end_highest == nullptr) {
    // Implicit H - RDKit uses the double bond partner as placeholder
    ranked_atoms[1] = begin_atom->getIdx();
  } else {
    ranked_atoms[1] = end_highest->getIdx();
  }

  // Compute and assign descriptor
  Descriptor desc = computeDoubleBondDescriptor(mol, bond, ranked_atoms);
  if (desc != Descriptor::NONE) {
    std::string desc_str = to_string(desc);
    if (!desc_str.empty()) {
      bond->setProp(common_properties::_CIPCode, desc_str);
    }
  }
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
