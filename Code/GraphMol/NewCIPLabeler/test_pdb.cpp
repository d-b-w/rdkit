//
//  Test NewCIPLabeler with PDB files
//
//  Usage: test_pdb <file.pdb>
//

#include <iostream>
#include <chrono>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "NewCIPLabeler.h"

using namespace RDKit;
using namespace std::chrono;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file.pdb>\n";
    return 1;
  }

  std::string pdb_file = argv[1];
  std::cout << "Reading PDB file: " << pdb_file << "\n";

  // Read PDB file
  auto start_parse = high_resolution_clock::now();
  std::unique_ptr<RWMol> mol(PDBFileToMol(pdb_file));
  auto end_parse = high_resolution_clock::now();

  if (!mol) {
    std::cerr << "Failed to parse PDB file\n";
    return 1;
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
    return 0;
  }

  // Wait for user to attach profiler
  std::cout << "\nReady to run CIP labeling.\n";
  std::cout << "Press Enter to continue (attach profiler now if needed)...";
  std::cin.get();

  // Run CIP labeling
  std::cout << "\nRunning CIP labeling...\n";
  auto start_cip = high_resolution_clock::now();

  try {
    NewCIPLabeler::assignCIPLabels(*mol);
  } catch (const NewCIPLabeler::MaxIterationsExceeded& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
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

  // Show first few labels as examples
  if (labeled_atoms > 0) {
    std::cout << "\nExample atom labels:\n";
    int shown = 0;
    for (auto& atom : mol->atoms()) {
      std::string code;
      if (atom->getPropIfPresent(common_properties::_CIPCode, code)) {
        std::cout << "  Atom " << atom->getIdx() << " ("
                  << atom->getSymbol() << "): " << code << "\n";
        if (++shown >= 10) break;
      }
    }
  }

  if (labeled_bonds > 0) {
    std::cout << "\nExample bond labels:\n";
    int shown = 0;
    for (auto& bond : mol->bonds()) {
      std::string code;
      if (bond->getPropIfPresent(common_properties::_CIPCode, code)) {
        std::cout << "  Bond " << bond->getIdx() << " ("
                  << bond->getBeginAtom()->getSymbol() << "-"
                  << bond->getEndAtom()->getSymbol() << "): " << code << "\n";
        if (++shown >= 10) break;
      }
    }
  }

  // Performance summary
  std::cout << "\nPerformance:\n";
  if (labeled_atoms + labeled_bonds > 0) {
    double us_per_center = cip_time.count() / static_cast<double>(labeled_atoms + labeled_bonds);
    std::cout << "  Average: " << us_per_center << " µs per stereocenter\n";
    std::cout << "  Rate: " << (1000000.0 / us_per_center) << " centers/sec\n";
  }

  return 0;
}
