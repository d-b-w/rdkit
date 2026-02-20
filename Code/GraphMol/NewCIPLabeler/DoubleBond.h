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
class Bond;

namespace NewCIPLabeler {

// Label a double bond with E/Z descriptor
void labelDoubleBond(ROMol& mol, Bond* bond, uint32_t max_iters);

}  // namespace NewCIPLabeler
}  // namespace RDKit
