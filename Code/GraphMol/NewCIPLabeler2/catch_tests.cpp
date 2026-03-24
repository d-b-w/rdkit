//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>

#include <catch2/catch_all.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>

#include "NewCIPLabeler.h"

using namespace RDKit;
using namespace RDKit::NewCIPLabeler2;

TEST_CASE("Stub - CIP computed flag set", "[newCIP2][stub]") {
  SECTION("Simple molecule") {
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    bool computed = false;
    mol->getPropIfPresent(common_properties::_CIPComputed, computed);
    CHECK(computed);
  }

  SECTION("Simple molecule with bitsets") {
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    boost::dynamic_bitset<> atoms(mol->getNumAtoms());
    atoms.set();
    boost::dynamic_bitset<> bonds(mol->getNumBonds());
    bonds.set();

    assignCIPLabels(*mol, atoms, bonds);

    bool computed = false;
    mol->getPropIfPresent(common_properties::_CIPComputed, computed);
    CHECK(computed);
  }
}

TEST_CASE("Stub - No labels assigned", "[newCIP2][stub]") {
  SECTION("Stub implementation doesn't assign labels") {
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto *center = mol->getAtomWithIdx(1);
    std::string cip_code;
    bool has_label =
        center->getPropIfPresent(common_properties::_CIPCode, cip_code);

    // Stub implementation doesn't assign labels
    CHECK_FALSE(has_label);
  }
}
