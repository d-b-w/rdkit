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

#include "Substituent.h"
#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

namespace NewCIPLabeler {

// Try to rank substituents using only constitutional rules (atomic number & mass)
// Returns true if all substituents are uniquely ranked
bool tryConstitutionalRanking(const ROMol& mol,
                              std::vector<Substituent>& subs);

// Apply full CIP ranking rules at a specific shell depth
// Returns true if all substituents are uniquely ranked
bool applyRankingRules(const ROMol& mol,
                       std::vector<Substituent>& subs,
                       uint32_t shell_depth);

// Check if the configuration represents pseudo-asymmetry
bool checkPseudoAsymmetry(const std::vector<Substituent>& subs);

}  // namespace NewCIPLabeler
}  // namespace RDKit
