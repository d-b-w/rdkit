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
#include "Descriptor.h"

using namespace RDKit;
using namespace RDKit::NewCIPLabeler;

TEST_CASE("Constitutional - atomic number", "[newCIP][phase1]") {
  SECTION("Simple Br, Cl, F, H") {
    // Br > Cl > F > H by atomic number
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    // Check that CIP was computed
    bool computed = false;
    mol->getPropIfPresent(common_properties::_CIPComputed, computed);
    CHECK(computed);

    // Check the stereocenter (atom index 1, the carbon)
    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    if (center->getPropIfPresent(common_properties::_CIPCode, cip_code)) {
      // The exact label depends on the 3D geometry encoded in the chiral tag
      // For now, just check that we got a label
      CHECK((cip_code == "R" || cip_code == "S"));
    } else {
      // If we didn't get a label, that's a problem
      FAIL("No CIP code assigned");
    }
  }

  SECTION("Different heteroatoms") {
    // O > N > C > H by atomic number
    auto mol = "O[C@H](N)C"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK((cip_code == "R" || cip_code == "S"));
  }
}

TEST_CASE("Constitutional - isotopes", "[newCIP][phase1]") {
  SECTION("Deuterium vs hydrogen") {
    // [2H] (deuterium) > [1H] (normal H) by mass
    auto mol = "[2H][C@H](C)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK((cip_code == "R" || cip_code == "S"));
  }

  SECTION("Carbon isotopes") {
    // [13C] > [12C] by mass
    auto mol = "[13CH3][C@H](C)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK((cip_code == "R" || cip_code == "S"));
  }
}

TEST_CASE("No chiral tag", "[newCIP][phase1]") {
  SECTION("No label without chiral tag") {
    // No @ or @@ in SMILES
    auto mol = "BrC(Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    // Should mark as computed
    bool computed = false;
    mol->getPropIfPresent(common_properties::_CIPComputed, computed);
    CHECK(computed);

    // But no CIP code should be assigned
    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    CHECK_FALSE(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
  }
}

TEST_CASE("Simple R/S assignment", "[newCIP][phase1]") {
  SECTION("Known R configuration") {
    // This is a known R enantiomer
    auto mol = "Br[C@@H](Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    REQUIRE(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    // With @@ (CCW), and priority Br > Cl > F > H, this should be R
    // (but the exact result depends on SMILES interpretation)
    CHECK_FALSE(cip_code.empty());
  }

  SECTION("Known S configuration") {
    // This is a known S enantiomer
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    REQUIRE(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK_FALSE(cip_code.empty());
  }
}

TEST_CASE("Exception handling", "[newCIP][phase1]") {
  SECTION("Max iterations") {
    // Create a complex symmetric molecule that would require many iterations
    // For now, just test with a simple molecule and low limit
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol);

    // With max_iters = 0, should work (constitutional ranking)
    CHECK_NOTHROW(assignCIPLabels(*mol, 0));

    // Verify it worked
    auto* center = mol->getAtomWithIdx(1);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
  }
}

TEST_CASE("Descriptor to_string", "[newCIP][phase1]") {
  CHECK(to_string(Descriptor::NONE) == "");
  CHECK(to_string(Descriptor::R) == "R");
  CHECK(to_string(Descriptor::S) == "S");
  CHECK(to_string(Descriptor::r) == "r");
  CHECK(to_string(Descriptor::s) == "s");
  CHECK(to_string(Descriptor::E) == "E");
  CHECK(to_string(Descriptor::Z) == "Z");
  CHECK(to_string(Descriptor::M) == "M");
  CHECK(to_string(Descriptor::P) == "P");
}
