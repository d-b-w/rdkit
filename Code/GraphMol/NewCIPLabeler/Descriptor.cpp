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

Descriptor descriptorFromString(const std::string& str) {
  if (str.empty()) return Descriptor::NONE;
  if (str == "R") return Descriptor::R;
  if (str == "S") return Descriptor::S;
  if (str == "r") return Descriptor::r;
  if (str == "s") return Descriptor::s;
  if (str == "E") return Descriptor::E;
  if (str == "Z") return Descriptor::Z;
  if (str == "M") return Descriptor::M;
  if (str == "P") return Descriptor::P;
  if (str == "m") return Descriptor::m;
  if (str == "p") return Descriptor::p;
  if (str == "ns") return Descriptor::ns;
  if (str == "seqCis") return Descriptor::seqCis;
  if (str == "seqTrans") return Descriptor::seqTrans;
  if (str == "SP-4") return Descriptor::SP_4;
  if (str == "TBPY-5") return Descriptor::TBPY_5;
  if (str == "OC-6") return Descriptor::OC_6;

  return Descriptor::UNKNOWN;
}

}  // namespace NewCIPLabeler
}  // namespace RDKit
