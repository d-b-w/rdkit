//
//  Test NewCIPLabeler with PDB files
//
//  Usage: test_pdb <file1.pdb> [file2.pdb ...]
//  Supports multiple files and glob patterns
//

//
// mmshare/test/testfiles/4axm.pdb
//
// /Users/dbn/builds/26-2/source/mmshare/test/testfiles/2j3n.pdb
// /Users/dbn/builds/26-2/source/mmshare/maestrolibs/src/structhierarchy/test/test_data/6m17.pdb
// /Users/dbn/builds/26-2/source/mmshare/test/testfiles/bioluminate/antibody/1baf.pdb
// /Users/dbn/builds/26-2/source/mmshare/test/testfiles/1qpe.pdb
// /Users/dbn/builds/26-2/source/mmshare/test/testfiles/ccd_mismatch.pdb
// /Users/dbn/builds/26-2/source/mmshare/test/schrodinger/structure/mmpdb/testsuite/1ida.pdb
// /Users/dbn/builds/26-2/source/mmshare/python/test/common_scripts/primex_polish_tests/2i35.pdb
// /Users/dbn/builds/26-2/source/mmshare/python/test/common_scripts/primex_polish_tests/test_data/2i35.pdb
// /Users/dbn/builds/26-2/source/mmshare/python/test/common_scripts/primex_polish_tests/test_data/PRIMEX-1174/3mbv.pdb
//


#include <iostream>
#include <chrono>
#include <string>
#include <exception>
#include <vector>
#include <algorithm>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include "NewCIPLabeler.h"

using namespace RDKit;
using namespace std::chrono;

// Track file performance and failures
struct FileResult {
  std::string filename;
  long long time_us;
  int stereocenters;
  bool success;
  std::string error;
};

std::vector<FileResult> g_results;

// using NewCIPLabeler::assignCIPLabels;
// using NewCIPLabeler::MaxIterationsExceeded;
using CIPLabeler::assignCIPLabels;
using CIPLabeler::MaxIterationsExceeded;

static bool processPDBFile(const std::string& pdb_file, FileResult& result) {
  result.filename = pdb_file;
  result.success = false;
  result.time_us = 0;
  result.stereocenters = 0;
  result.error = "";

  std::cout << "\n" << std::string(80, '=') << "\n";
  std::cout << "Reading PDB file: " << pdb_file << "\n";

  // Read PDB file
  auto start_parse = high_resolution_clock::now();
  std::unique_ptr<RWMol> mol;

  try {
    mol.reset(PDBFileToMol(pdb_file));
  } catch (const std::exception& e) {
    result.error = std::string("Parse error: ") + e.what();
    std::cerr << "ERROR: Failed to parse PDB file: " << e.what() << "\n";
    std::cerr << "Skipping this file.\n";
    return false;
  }

  auto end_parse = high_resolution_clock::now();

  if (!mol) {
    result.error = "Parse error: returned null";
    std::cerr << "ERROR: Failed to parse PDB file (returned null)\n";
    std::cerr << "Skipping this file.\n";
    return false;
  }

  auto parse_time = duration_cast<milliseconds>(end_parse - start_parse);
  std::cout << "Parsed in " << parse_time.count() << " ms\n";

  // Molecule statistics
  // std::cout << "\nMolecule statistics:\n";
  // std::cout << "  Atoms: " << mol->getNumAtoms() << "\n";
  // std::cout << "  Bonds: " << mol->getNumBonds() << "\n";

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

  // std::cout << "  Tetrahedral centers: " << tetrahedral_count << "\n";
  // std::cout << "  Stereo double bonds: " << double_bond_count << "\n";

  result.stereocenters = tetrahedral_count + double_bond_count;

  if (tetrahedral_count == 0 && double_bond_count == 0) {
    std::cout << "\nNo stereocenters found - nothing to label\n";
    result.success = true;
    return true;  // Still successful, just nothing to do
  }

  // Run CIP labeling
  std::cout << "\nRunning CIP labeling...\n";
  auto start_cip = high_resolution_clock::now();

  try {
    assignCIPLabels(*mol);
  } catch (const MaxIterationsExceeded& e) {
    auto end_cip = high_resolution_clock::now();
    auto cip_time = duration_cast<microseconds>(end_cip - start_cip);
    result.time_us = cip_time.count();
    result.error = std::string("MaxIterationsExceeded: ") + e.what();
    std::cerr << "ERROR: " << e.what() << "\n";
    std::cerr << " after " << cip_time.count() << " µs (" << (cip_time.count() / 1000.0) << " ms)\n";
    return false;
  } catch (const std::exception& e) {
    auto end_cip = high_resolution_clock::now();
    auto cip_time = duration_cast<microseconds>(end_cip - start_cip);
    result.time_us = cip_time.count();
    result.error = std::string("Exception: ") + e.what();
    std::cerr << "ERROR: " << e.what() << "\n";
    std::cerr << " after " << cip_time.count() << " µs (" << (cip_time.count() / 1000.0) << " ms)\n";
    return false;
  }

  auto end_cip = high_resolution_clock::now();
  auto cip_time = duration_cast<microseconds>(end_cip - start_cip);

  result.time_us = cip_time.count();
  result.success = true;

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

  // std::cout << "\nResults:\n";
  // std::cout << "  Labeled atoms: " << labeled_atoms << " / " << tetrahedral_count << "\n";
  // std::cout << "  Labeled bonds: " << labeled_bonds << " / " << double_bond_count << "\n";

  // Show atom labels
  if (labeled_atoms > 0 && false) {
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

  if (labeled_bonds > 0 && false) {
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
  // std::cout << "\nPerformance:\n";
  if (labeled_atoms + labeled_bonds > 0) {
    double us_per_center = cip_time.count() / static_cast<double>(labeled_atoms + labeled_bonds);
    std::cout << "  Average: " << us_per_center << " µs per stereocenter\n";
    // bold if slow
    if (cip_time.count() > 100000) {
        std::cout << "\e[1m";
    }
    std::cout << "  Total: " << cip_time.count() << " µs";
    if (cip_time.count() > 100000) {
        std::cout << "\e[0m";
    }
    // std::cout << "\n  Rate: " << (1000000.0 / us_per_center) << " centers/sec\n";
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
  // std::cin.get();

  for (int i = 1; i < argc; ++i) {
    FileResult result;
    try {
      auto success = processPDBFile(argv[i], result);
      g_results.push_back(result);
      if (success) {
        successful++;
      } else {
        failed++;
      }
    } catch (const std::exception& e) {
      result.filename = argv[i];
      result.success = false;
      result.error = std::string("FATAL: ") + e.what();
      g_results.push_back(result);
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

  // Show failed files (CIP labeling errors only, not parse errors)
  std::vector<FileResult> cip_failures;
  for (const auto& r : g_results) {
    if (!r.success && r.error.find("Parse error:") != 0) {
      cip_failures.push_back(r);
    }
  }

  if (!cip_failures.empty()) {
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << "CIP labeling failures:\n";
    for (const auto& r : cip_failures) {
      std::cout << "  " << r.filename << "\n";
      std::cout << "    Error: " << r.error << "\n";
      if (r.time_us > 0) {
        std::cout << "    Time before failure: " << r.time_us << " µs ("
                  << (r.time_us / 1000.0) << " ms)\n";
      }
    }
  }

  // Show slowest files (top 10)
  if (successful > 0) {
    std::vector<FileResult> successful_results;
    for (const auto& r : g_results) {
      if (r.success && r.stereocenters > 0) {
        successful_results.push_back(r);
      }
    }

    if (!successful_results.empty()) {
      std::sort(successful_results.begin(), successful_results.end(),
                [](const FileResult& a, const FileResult& b) {
                  return a.time_us > b.time_us;
                });

      std::cout << "\n" << std::string(80, '-') << "\n";
      std::cout << "Slowest files (top " << std::min(10, (int)successful_results.size()) << "):\n";
      int shown = 0;
      for (const auto& r : successful_results) {
        if (shown >= 10) break;
        double ms = r.time_us / 1000.0;
        double us_per_center = r.time_us / static_cast<double>(r.stereocenters);
        std::cout << "  " << (shown + 1) << ". " << r.filename << "\n";
        std::cout << "     " << r.time_us << " µs (" << ms << " ms) for "
                  << r.stereocenters << " centers (" << us_per_center
                  << " µs/center)\n";
        shown++;
      }
    }
  }

  return (failed == total_files) ? 1 : 0;
}
