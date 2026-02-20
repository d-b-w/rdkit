//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Descriptor.h"
#include <stdexcept>

namespace RDKit {
namespace NewCIPLabeler {

std::string to_string(Descriptor desc) {
  switch (desc) {
    case Descriptor::NONE:
      return "";
    case Descriptor::UNKNOWN:
      return "";
    case Descriptor::ns:
      return "ns";
    case Descriptor::R:
      return "R";
    case Descriptor::S:
      return "S";
    case Descriptor::r:
      return "r";
    case Descriptor::s:
      return "s";
    case Descriptor::E:
      return "E";
    case Descriptor::Z:
      return "Z";
    case Descriptor::seqCis:
      return "seqCis";
    case Descriptor::seqTrans:
      return "seqTrans";
    case Descriptor::M:
      return "M";
    case Descriptor::P:
      return "P";
    case Descriptor::m:
      return "m";
    case Descriptor::p:
      return "p";
    case Descriptor::SP_4:
      return "SP-4";
    case Descriptor::TBPY_5:
      return "TBPY-5";
    case Descriptor::OC_6:
      return "OC-6";
    default:
      throw std::runtime_error("Unknown descriptor");
  }
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
