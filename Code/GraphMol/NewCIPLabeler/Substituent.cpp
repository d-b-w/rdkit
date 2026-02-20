//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Substituent.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>

namespace RDKit {
namespace NewCIPLabeler {

void Substituent::expandNextShell(const ROMol& mol,
                                  boost::dynamic_bitset<>& visited) {
  // Handle implicit hydrogen (root_atom is nullptr)
  if (root_atom == nullptr) {
    if (shells.empty()) {
      // Shell 0 for implicit H: single hydrogen atom (conceptual)
      AtomShell shell{0, {nullptr}, {1}};
      shells.push_back(std::move(shell));
    }
    // Hydrogen has no further shells
    return;
  }

  uint32_t new_shell_idx = static_cast<uint32_t>(shells.size());
  AtomShell new_shell{new_shell_idx, {}, {}};

  if (new_shell_idx == 0) {
    // Shell 0: just the root atom
    new_shell.atoms.push_back(root_atom);
    new_shell.multiplicities.push_back(1);
    visited.set(root_atom->getIdx());
  } else {
    // Expand from previous shell
    const auto& prev = shells[new_shell_idx - 1];

    for (size_t i = 0; i < prev.atoms.size(); ++i) {
      const Atom* atom = prev.atoms[i];
      if (atom == nullptr) {
        continue;  // Skip nullptr (shouldn't happen except for H)
      }

      uint32_t parent_mult = prev.multiplicities[i];

      for (const auto& bond : mol.atomBonds(atom)) {
        const Atom* nbr = bond->getOtherAtom(atom);

        if (visited[nbr->getIdx()]) {
          // Ring closure - add as duplicate node with multiplicity 0
          new_shell.atoms.push_back(nbr);
          new_shell.multiplicities.push_back(0);
          continue;
        }

        // New atom - add with bond order multiplicity
        int bond_order = 1;
        switch (bond->getBondType()) {
          case Bond::SINGLE:
            bond_order = 1;
            break;
          case Bond::DOUBLE:
            bond_order = 2;
            break;
          case Bond::TRIPLE:
            bond_order = 3;
            break;
          case Bond::QUADRUPLE:
            bond_order = 4;
            break;
          case Bond::AROMATIC:
            bond_order = 1;  // Treat aromatic as 1 for CIP
            break;
          case Bond::DATIVE:
            bond_order = 0;  // Dative bonds have 0 order for CIP
            break;
          default:
            bond_order = 1;
        }

        new_shell.atoms.push_back(nbr);
        new_shell.multiplicities.push_back(parent_mult * bond_order);
        visited.set(nbr->getIdx());
      }
    }
  }

  shells.push_back(std::move(new_shell));
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
