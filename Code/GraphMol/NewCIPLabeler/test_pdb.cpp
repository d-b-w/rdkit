//
//  Test NewCIPLabeler with PDB files
//
//  Usage: test_pdb <file1.pdb> [file2.pdb ...]
//  Supports multiple files and glob patterns
//

#include <iostream>
#include <chrono>
#include <string>
#include <exception>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include "NewCIPLabeler.h"

using namespace RDKit;
using namespace std::chrono;

using NewCIPLabeler::assignCIPLabels;
using NewCIPLabeler::MaxIterationsExceeded;
// using CIPLabeler::assignCIPLabels;
// using CIPLabeler::MaxIterationsExceeded;

bool processPDBFile(const std::string& pdb_file) {
  std::cout << "\n" << std::string(80, '=') << "\n";
  std::cout << "Reading PDB file: " << pdb_file << "\n";

  // Read PDB file
  auto start_parse = high_resolution_clock::now();
  std::unique_ptr<RWMol> mol;

  try {
    mol.reset(PDBFileToMol(pdb_file));
  } catch (const std::exception& e) {
    std::cerr << "ERROR: Failed to parse PDB file: " << e.what() << "\n";
    std::cerr << "Skipping this file.\n";
    return false;
  }

  auto end_parse = high_resolution_clock::now();

  if (!mol) {
    std::cerr << "ERROR: Failed to parse PDB file (returned null)\n";
    std::cerr << "Skipping this file.\n";
    return false;
  }

  auto parse_time = duration_cast<milliseconds>(end_parse - start_parse);
  std::cout << "Parsed in " << parse_time.count() << " ms\n";

  // Molecule statistics
  std::cout << "\nMolecule statistics:\n";
  std::cout << "  Atoms: " << mol->getNumAtoms() << "\n";
  std::cout << "  Bonds: " << mol->getNumBonds() << "\n";

  // Count potential stereocenters
  int tetrahedral_count = 0;
  for (auto& atom : mol->atoms()) {
    auto tag = atom->getChiralTag();
    if (tag == Atom::CHI_TETRAHEDRAL_CW || tag == Atom::CHI_TETRAHEDRAL_CCW) {
      tetrahedral_count++;
    }
  }

  int double_bond_count = 0;
  for (auto& bond : mol->bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      auto stereo = bond->getStereo();
      if (stereo == Bond::STEREOCIS || stereo == Bond::STEREOTRANS ||
          stereo == Bond::STEREOE || stereo == Bond::STEREOZ) {
        double_bond_count++;
      }
    }
  }

  std::cout << "  Tetrahedral centers: " << tetrahedral_count << "\n";
  std::cout << "  Stereo double bonds: " << double_bond_count << "\n";

  if (tetrahedral_count == 0 && double_bond_count == 0) {
    std::cout << "\nNo stereocenters found - nothing to label\n";
    return true;  // Still successful, just nothing to do
  }

  // Run CIP labeling
  std::cout << "\nRunning CIP labeling...\n";
  auto start_cip = high_resolution_clock::now();

  try {
    assignCIPLabels(*mol);
  } catch (const MaxIterationsExceeded& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return false;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return false;
  }

  auto end_cip = high_resolution_clock::now();
  auto cip_time = duration_cast<microseconds>(end_cip - start_cip);

  std::cout << "CIP labeling completed in " << cip_time.count() << " µs ("
            << (cip_time.count() / 1000.0) << " ms)\n";

  // Count labeled centers
  int labeled_atoms = 0;
  for (auto& atom : mol->atoms()) {
    std::string code;
    if (atom->getPropIfPresent(common_properties::_CIPCode, code)) {
      labeled_atoms++;
    }
  }

  int labeled_bonds = 0;
  for (auto& bond : mol->bonds()) {
    std::string code;
    if (bond->getPropIfPresent(common_properties::_CIPCode, code)) {
      labeled_bonds++;
    }
  }

  std::cout << "\nResults:\n";
  std::cout << "  Labeled atoms: " << labeled_atoms << " / " << tetrahedral_count << "\n";
  std::cout << "  Labeled bonds: " << labeled_bonds << " / " << double_bond_count << "\n";

  // Show atom labels
  if (labeled_atoms > 0) {
    std::cout << "\nAtom labels";
    if (labeled_atoms > 5) {
      std::cout << " (showing first 5 of " << labeled_atoms << ")";
    }
    std::cout << ":\n";

    int shown = 0;
    for (auto& atom : mol->atoms()) {
      std::string code;
      if (atom->getPropIfPresent(common_properties::_CIPCode, code)) {
        std::cout << "  Atom " << atom->getIdx() << " ("
                  << atom->getSymbol() << "): " << code << "\n";
        if (++shown >= 5) break;
      }
    }

    if (labeled_atoms > 5) {
      std::cout << "  ... and " << (labeled_atoms - 5) << " more\n";
    }
  }

  if (labeled_bonds > 0) {
    std::cout << "\nBond labels";
    if (labeled_bonds > 5) {
      std::cout << " (showing first 5 of " << labeled_bonds << ")";
    }
    std::cout << ":\n";

    int shown = 0;
    for (auto& bond : mol->bonds()) {
      std::string code;
      if (bond->getPropIfPresent(common_properties::_CIPCode, code)) {
        std::cout << "  Bond " << bond->getIdx() << " ("
                  << bond->getBeginAtom()->getSymbol() << "-"
                  << bond->getEndAtom()->getSymbol() << "): " << code << "\n";
        if (++shown >= 5) break;
      }
    }

    if (labeled_bonds > 5) {
      std::cout << "  ... and " << (labeled_bonds - 5) << " more\n";
    }
  }

  // Performance summary
  std::cout << "\nPerformance:\n";
  if (labeled_atoms + labeled_bonds > 0) {
    double us_per_center = cip_time.count() / static_cast<double>(labeled_atoms + labeled_bonds);
    std::cout << "  Average: " << us_per_center << " µs per stereocenter\n";
    std::cout << "  Total: " << cip_time.count() << " µs\n";
    std::cout << "  Rate: " << (1000000.0 / us_per_center) << " centers/sec\n";
  }

  return true;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file1.pdb> [file2.pdb ...]\n";
    std::cerr << "Process one or more PDB files (use glob patterns)\n";
    return 1;
  }

  int total_files = argc - 1;
  int successful = 0;
  int failed = 0;

  std::cout << "Processing " << total_files << " PDB file(s)...\n";
  // Wait for user to attach profiler (only for first file)
  std::cout << "\nReady to run CIP labeling.\n";
  std::cout << "Press Enter to continue (attach profiler now if needed)...";
  std::cin.get();

  for (int i = 1; i < argc; ++i) {
    try {
      bool success = processPDBFile(argv[i]);
      if (success) {
        successful++;
      } else {
        failed++;
      }
    } catch (const std::exception& e) {
      std::cerr << "\nFATAL ERROR processing " << argv[i] << ": " << e.what() << "\n";
      std::cerr << "Skipping this file.\n";
      failed++;
    }
  }

  // Summary
  std::cout << "\n" << std::string(80, '=') << "\n";
  std::cout << "Summary:\n";
  std::cout << "  Total files: " << total_files << "\n";
  std::cout << "  Successful: " << successful << "\n";
  std::cout << "  Failed: " << failed << "\n";

  return (failed == total_files) ? 1 : 0;
}
