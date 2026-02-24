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

#include <RDGeneral/export.h>
#include <cstdint>

namespace RDKit {

class ROMol;
class Atom;

namespace NewCIPLabeler {

// Label a tetrahedral stereocenter with R/S descriptor
// use_rule5: enable CIP Rule 5 for pseudo-asymmetry detection
void labelTetrahedralCenter(ROMol& mol, Atom* center, uint32_t max_iters, bool use_rule5 = true);

}  // namespace NewCIPLabeler
}  // namespace RDKit
