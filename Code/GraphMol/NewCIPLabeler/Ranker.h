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
class Atom;

namespace NewCIPLabeler {

// Core ranking function - attempts to uniquely rank all substituents
// using lazy shell expansion up to max_shells depth
CenterRanking rankSubstituents(const ROMol& mol,
                               const Atom* center,
                               std::vector<Substituent>& subs,
                               uint32_t max_shells);

// Build final ranking result from ranked substituents
CenterRanking buildRankingResult(const std::vector<Substituent>& subs,
                                 bool is_pseudo);

}  // namespace NewCIPLabeler
}  // namespace RDKit
