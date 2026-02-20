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

#include <string>
#include <RDGeneral/export.h>

namespace RDKit {
namespace NewCIPLabeler {

// CIP descriptor types following the algorithmic paper
enum class Descriptor {
  NONE,      // Unspecified
  UNKNOWN,   // Cannot determine
  ns,        // Unspecified other

  // Tetrahedral
  R, S,      // Chiral
  r, s,      // Pseudo-chiral

  // Double bond
  E, Z,      // Entgegen/Zusammen
  seqCis, seqTrans,  // Sequential

  // Atropisomer
  M, P,      // Axial chirality
  m, p,      // Pseudo-axial

  // Coordination (future)
  SP_4, TBPY_5, OC_6
};

// Convert descriptor to string for storage in molecule properties
RDKIT_NEWCIPLABELER_EXPORT std::string to_string(Descriptor desc);

}  // namespace NewCIPLabeler
}  // namespace RDKit
