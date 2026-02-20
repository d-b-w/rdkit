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

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/export.h>
#include <stdexcept>

namespace RDKit {

class ROMol;

namespace NewCIPLabeler {

/**
 * Exception thrown when maximum iterations exceeded during CIP label calculation.
 * Very symmetrical molecules can cause excessive processing (e.g. dodecahedrane).
 */
class RDKIT_NEWCIPLABELER_EXPORT MaxIterationsExceeded
    : public std::runtime_error {
 public:
  explicit MaxIterationsExceeded()
      : std::runtime_error("Max Iterations Exceeded in CIP label calculation") {}
};

/**
 * Calculate stereochemical labels using high-performance CIP implementation.
 *
 * This is a performance-optimized implementation that uses lazy shell expansion
 * to dramatically improve speed over the standard CIPLabeler. The algorithm
 * follows the paper:
 *
 * Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A.,
 * Redkin, D. Algorithmic Analysis of Cahn--Ingold--Prelog Rules of
 * Stereochemistry: Proposals for Revised Rules and a Guide for Machine
 * Implementation. J. Chem. Inf. Model. 2018, 58, 1755-1765.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \param maxRecursiveIterations - maximum shell expansions (0 = unlimited)
 *
 *   \note only atoms with chiral tags and double bonds with proper
 *          bond directions will be labelled.
 *   \note Labels will be stored under the common_properties::_CIPCode
 *          property of the relevant atoms/bonds.
 */
RDKIT_NEWCIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol, unsigned int maxRecursiveIterations = 0);

/**
 * Overload that allows selecting which atoms and/or bonds will be labeled.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \param atoms - bitset with the atom indexes to be labeled.
 *
 *   \param bonds - bitset with the bond indexes to be labeled.
 *
 *   \param maxRecursiveIterations - maximum shell expansions (0 = unlimited)
 */
RDKIT_NEWCIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol, const boost::dynamic_bitset<> &atoms,
    const boost::dynamic_bitset<> &bonds,
    unsigned int maxRecursiveIterations = 0);

}  // namespace NewCIPLabeler
}  // namespace RDKit
