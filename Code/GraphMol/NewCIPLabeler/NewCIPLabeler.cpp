//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "NewCIPLabeler.h"
#include "Tetrahedral.h"
#include "DoubleBond.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <vector>
#include <algorithm>

namespace RDKit {
namespace NewCIPLabeler {

namespace {

// Find atoms that might be stereocenters (have chiral tags)
std::vector<Atom*> findStereoAtoms(ROMol& mol,
                                   const boost::dynamic_bitset<>& filter) {
  std::vector<Atom*> result;

  for (auto& atom : mol.atoms()) {
    unsigned int idx = atom->getIdx();
    if (!filter[idx]) {
      continue;
    }

    auto tag = atom->getChiralTag();
    if (tag == Atom::CHI_TETRAHEDRAL_CW || tag == Atom::CHI_TETRAHEDRAL_CCW) {
      result.push_back(atom);
    }
  }

  return result;
}

// Find bonds that might have stereochemistry (double bonds with directional info)
std::vector<Bond*> findStereoBonds(ROMol& mol,
                                   const boost::dynamic_bitset<>& filter) {
  std::vector<Bond*> result;

  for (auto& bond : mol.bonds()) {
    unsigned int idx = bond->getIdx();
    if (!filter[idx]) {
      continue;
    }

    // Check for double bond with stereo info
    if (bond->getBondType() == Bond::DOUBLE) {
      auto stereo = bond->getStereo();
      if (stereo == Bond::STEREOCIS || stereo == Bond::STEREOTRANS ||
          stereo == Bond::STEREOE || stereo == Bond::STEREOZ) {
        result.push_back(bond);
      }
    }
  }

  return result;
}

// Estimate complexity of a stereocenter for sorting
uint32_t estimateComplexity(const ROMol& mol, const Atom* atom) {
  uint32_t complexity = 0;

  // Ring membership increases complexity
  const auto& ringInfo = mol.getRingInfo();
  complexity += ringInfo->numAtomRings(atom->getIdx()) * 10;

  // Heteroatoms in neighborhood increase complexity
  for (const auto& nbr : mol.atomNeighbors(atom)) {
    if (nbr->getAtomicNum() > 6) {
      complexity += 5;
    }
  }

  // Count neighbor atomic numbers to detect symmetry
  std::map<int, int> z_counts;
  for (const auto& nbr : mol.atomNeighbors(atom)) {
    z_counts[nbr->getAtomicNum()]++;
  }

  // Multiple same-Z neighbors indicate potential symmetry
  for (const auto& [z, count] : z_counts) {
    if (count > 1) {
      complexity += count * 3;
    }
  }

  return complexity;
}

// Sort centers by estimated complexity (simple first)
std::vector<Atom*> sortByComplexity(const std::vector<Atom*>& atoms,
                                    const ROMol& mol) {
  std::vector<std::pair<uint32_t, Atom*>> scored;
  scored.reserve(atoms.size());

  for (auto* atom : atoms) {
    uint32_t score = estimateComplexity(mol, atom);
    scored.emplace_back(score, atom);
  }

  // Sort by score (ascending = simple first)
  std::sort(scored.begin(), scored.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

  std::vector<Atom*> result;
  result.reserve(atoms.size());
  for (const auto& [score, atom] : scored) {
    result.push_back(atom);
  }

  return result;
}

}  // namespace

void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations) {
  boost::dynamic_bitset<> all_atoms(mol.getNumAtoms());
  boost::dynamic_bitset<> all_bonds(mol.getNumBonds());
  all_atoms.set();
  all_bonds.set();

  assignCIPLabels(mol, all_atoms, all_bonds, maxRecursiveIterations);
}

void assignCIPLabels(ROMol &mol,
                     const boost::dynamic_bitset<> &atoms,
                     const boost::dynamic_bitset<> &bonds,
                     unsigned int maxRecursiveIterations) {
  // Collect potential stereocenters
  auto stereo_atoms = findStereoAtoms(mol, atoms);
  auto stereo_bonds = findStereoBonds(mol, bonds);

  // Sort by complexity (simple first)
  stereo_atoms = sortByComplexity(stereo_atoms, mol);

  // Process each stereocenter
  for (auto* atom : stereo_atoms) {
    try {
      labelTetrahedralCenter(mol, atom, maxRecursiveIterations);
    } catch (const MaxIterationsExceeded&) {
      // Let it propagate up
      throw;
    } catch (const std::exception& e) {
      // Silently skip atoms that fail labeling
      // This matches old CIPLabeler behavior
      continue;
    }
  }

  // Process double bonds
  for (auto* bond : stereo_bonds) {
    try {
      labelDoubleBond(mol, bond, maxRecursiveIterations);
    } catch (const MaxIterationsExceeded&) {
      // Let it propagate up
      throw;
    } catch (const std::exception& e) {
      // Silently skip bonds that fail labeling
      continue;
    }
  }

  // Mark molecule as having CIP computed
  mol.setProp(common_properties::_CIPComputed, true);
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
