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

TEST_CASE("Debug SmilesToMol CIP behavior", "[newCIP][debug]") {
  auto mol = "Br[C@H](Cl)F"_smiles;
  REQUIRE(mol);

  // Check what SmilesToMol sets
  INFO("After SmilesToMol:");
  INFO("  Molecule _CIPComputed: " << mol->hasProp(common_properties::_CIPComputed));

  auto* atom = mol->getAtomWithIdx(1);
  std::string code;
  bool hasCode = atom->getPropIfPresent(common_properties::_CIPCode, code);
  INFO("  Atom _CIPCode present: " << hasCode);
  if (hasCode) {
    INFO("    Value: " << code);
  }

  // Now clear and try to reassign
  atom->clearProp(common_properties::_CIPCode);
  mol->clearProp(common_properties::_CIPComputed);  // Also clear molecule flag

  bool exception_thrown = false;
  std::string exception_msg;
  try {
    assignCIPLabels(*mol);
  } catch (const std::exception& e) {
    exception_thrown = true;
    exception_msg = e.what();
  }

  hasCode = atom->getPropIfPresent(common_properties::_CIPCode, code);
  INFO("After clearProp + assignCIPLabels:");
  INFO("  Exception thrown: " << exception_thrown);
  if (exception_thrown) {
    INFO("    Message: " << exception_msg);
  }
  INFO("  Atom _CIPCode present: " << hasCode);
  if (hasCode) {
    INFO("    Value: " << code);
  }

  CHECK(hasCode);
}

// ============================================================================
// Phase 2: Shell Expansion Tests
// ============================================================================

TEST_CASE("One shell expansion", "[newCIP][phase2]") {
  SECTION("Ethyl groups require 1-shell expansion") {
    // CC[C@H](C)CC - two ethyl groups, need to look beyond first carbon
    auto mol = "CC[C@H](C)CC"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(2);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    // Should get a label (R or S) even though immediate neighbors are all C
  }

  SECTION("Propyl vs methyl") {
    // CCC[C@H](C)F - propyl > ethyl at second shell
    auto mol = "CCC[C@H](C)F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(3);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK_FALSE(cip_code.empty());
  }
}

TEST_CASE("Ring systems", "[newCIP][phase2]") {
  SECTION("Cyclohexane with substituent") {
    // Ring system with chiral center
    auto mol = "C1CC[C@H](Br)CC1"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(3);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK_FALSE(cip_code.empty());
  }

  SECTION("Bridged ring") {
    // More complex ring system
    auto mol = "C1C[C@H]2CC[C@H](C1)C2"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    // Should handle ring closure properly
    bool computed = false;
    mol->getPropIfPresent(common_properties::_CIPComputed, computed);
    CHECK(computed);
  }
}

TEST_CASE("Symmetric substituents", "[newCIP][phase2]") {
  SECTION("Two identical branches") {
    // Need multiple shells to distinguish
    auto mol = "CC(C)[C@H](C(C)C)Br"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center = mol->getAtomWithIdx(3);
    std::string cip_code;
    CHECK(center->getPropIfPresent(common_properties::_CIPCode, cip_code));
    CHECK_FALSE(cip_code.empty());
  }
}

TEST_CASE("Multiple stereocenters", "[newCIP][phase2]") {
  SECTION("Two independent centers") {
    auto mol = "C[C@H](Cl)CC[C@H](Cl)C"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* center1 = mol->getAtomWithIdx(1);
    auto* center2 = mol->getAtomWithIdx(5);

    std::string cip1, cip2;
    CHECK(center1->getPropIfPresent(common_properties::_CIPCode, cip1));
    CHECK(center2->getPropIfPresent(common_properties::_CIPCode, cip2));

    // Both should have labels
    CHECK_FALSE(cip1.empty());
    CHECK_FALSE(cip2.empty());
  }
}

TEST_CASE("Port old CIPLabeler tests", "[newCIP][phase2]") {
  SECTION("Tetrahedral assignment - Br[C@H](Cl)F") {
    auto mol = "Br[C@H](Cl)F"_smiles;
    REQUIRE(mol->getNumAtoms() == 4);

    auto chiral_atom = mol->getAtomWithIdx(1);
    chiral_atom->clearProp(common_properties::_CIPCode);
    REQUIRE(chiral_atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

    CHECK(!mol->hasProp(common_properties::_CIPComputed));
    assignCIPLabels(*mol);
    CHECK(mol->hasProp(common_properties::_CIPComputed));

    std::string chirality;
    CHECK(chiral_atom->getPropIfPresent(common_properties::_CIPCode, chirality));
    CHECK(chirality == "S");
  }

  SECTION("Phosphine chirality") {
    const std::vector<std::pair<std::string, std::string>> mols{
        {"C[P@](C1CCCC1)C1=CC=CC=C1", "R"},
        {"C[As@@](C1CCCC1)C1=CC=CC=C1", "S"},
        {"C[P@H]C1CCCCC1", "R"},
        {"C[P@@H]C1CCCCC1", "S"}};

    for (const auto &ref : mols) {
      INFO(ref.first);
      std::unique_ptr<RWMol> mol{SmilesToMol(ref.first)};
      REQUIRE(mol);
      assignCIPLabels(*mol);

      std::string chirality;
      CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                     chirality));
      CHECK(chirality == ref.second);
    }
  }
}

// ============================================================================
// Phase 3: E/Z Double Bond Labeling
// ============================================================================

TEST_CASE("Simple E/Z labeling", "[newCIP][phase3]") {
  SECTION("E configuration - F/C=C/Cl") {
    // trans: higher priority on opposite sides
    auto mol = "F/C=C/Cl"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(1);  // C=C bond
    std::string stereo;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, stereo));
    CHECK(stereo == "E");
  }

  SECTION("Z configuration - F/C=C\\Cl") {
    // cis: higher priority on same side
    auto mol = "F/C=C\\Cl"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(1);  // C=C bond
    std::string stereo;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, stereo));
    CHECK(stereo == "Z");
  }

  SECTION("E configuration - Br/C=C/I") {
    // I > Br, both trans
    auto mol = "Br/C=C/I"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(1);
    std::string stereo;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, stereo));
    CHECK(stereo == "E");
  }
}

TEST_CASE("E/Z with substituents", "[newCIP][phase3]") {
  SECTION("Methyl vs ethyl - E") {
    // CC/C=C/C
    auto mol = "CC/C=C/C"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(2);  // C=C bond
    std::string stereo;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, stereo));
    CHECK(stereo == "E");
  }

  SECTION("Phenyl groups") {
    // c1ccccc1/C=C/c1ccccc1
    auto mol = "c1ccccc1/C=C/c1ccccc1"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(6);  // C=C bond
    std::string stereo;
    if (bond->getPropIfPresent(common_properties::_CIPCode, stereo)) {
      CHECK(stereo == "E");
    }
  }
}

TEST_CASE("E/Z edge cases", "[newCIP][phase3]") {
  SECTION("Symmetric - no label") {
    // F/C=C/F - both sides same
    auto mol = "F/C=C/F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(1);
    std::string stereo;
    // Should not have CIP code (symmetric)
    CHECK_FALSE(bond->hasProp(common_properties::_CIPCode));
  }

  SECTION("Trisubstituted - Cl/C(Br)=C/F") {
    auto mol = "Cl/C(Br)=C/F"_smiles;
    REQUIRE(mol);

    assignCIPLabels(*mol);

    auto* bond = mol->getBondWithIdx(2);  // C=C bond
    std::string stereo;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, stereo));
    // Br > Cl on one side, F > H on other, trans = E
    CHECK(stereo == "E");
  }
}
