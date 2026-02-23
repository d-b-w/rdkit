// Quick benchmark to check for obvious performance issues
// Compile: g++ -std=c++17 -O2 -I... quick_benchmark.cpp -lRDKit...

#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/RDKitBase.h>
#include "NewCIPLabeler.h"

using namespace RDKit;
using namespace std::chrono;

struct BenchmarkCase {
  std::string name;
  std::string smiles;
  int iterations;
};

void benchmark(const std::string& name, const std::string& smiles, int iterations) {
  std::vector<std::unique_ptr<ROMol>> mols;

  // Pre-parse molecules to isolate labeling time
  for (int i = 0; i < iterations; ++i) {
    mols.push_back(std::unique_ptr<ROMol>(SmilesToMol(smiles)));
    if (!mols.back()) {
      std::cerr << "Failed to parse: " << smiles << "\n";
      return;
    }
  }

  auto start = high_resolution_clock::now();

  for (auto& mol : mols) {
    NewCIPLabeler::assignCIPLabels(*mol);
  }

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);

  double avg_us = duration.count() / static_cast<double>(iterations);

  std::cout << name << ":\n";
  std::cout << "  SMILES: " << smiles << "\n";
  std::cout << "  Total: " << duration.count() << " µs for " << iterations << " iterations\n";
  std::cout << "  Average: " << avg_us << " µs per molecule\n";
  std::cout << "  Rate: " << (1000000.0 / avg_us) << " molecules/sec\n\n";
}

void benchmarkPDB(const std::string& pdb_file, int iterations) {
  std::cout << "PDB File: " << pdb_file << "\n";

  // Read PDB file once
  std::unique_ptr<RWMol> template_mol(PDBFileToMol(pdb_file));
  if (!template_mol) {
    std::cerr << "  Failed to parse PDB file\n\n";
    return;
  }

  std::cout << "  Atoms: " << template_mol->getNumAtoms()
            << ", Bonds: " << template_mol->getNumBonds() << "\n";

  // Count stereocenters
  int tetrahedral = 0, double_bonds = 0;
  for (auto& atom : template_mol->atoms()) {
    auto tag = atom->getChiralTag();
    if (tag == Atom::CHI_TETRAHEDRAL_CW || tag == Atom::CHI_TETRAHEDRAL_CCW) {
      tetrahedral++;
    }
  }
  for (auto& bond : template_mol->bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      auto stereo = bond->getStereo();
      if (stereo == Bond::STEREOCIS || stereo == Bond::STEREOTRANS ||
          stereo == Bond::STEREOE || stereo == Bond::STEREOZ) {
        double_bonds++;
      }
    }
  }

  std::cout << "  Stereocenters: " << tetrahedral << " tetrahedral, "
            << double_bonds << " double bonds\n";

  if (tetrahedral + double_bonds == 0) {
    std::cout << "  No stereocenters - skipping\n\n";
    return;
  }

  // Create copies for benchmarking
  std::vector<std::unique_ptr<RWMol>> mols;
  for (int i = 0; i < iterations; ++i) {
    mols.push_back(std::unique_ptr<RWMol>(new RWMol(*template_mol)));
  }

  auto start = high_resolution_clock::now();

  for (auto& mol : mols) {
    NewCIPLabeler::assignCIPLabels(*mol);
  }

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);

  double avg_us = duration.count() / static_cast<double>(iterations);

  std::cout << "  Total: " << duration.count() << " µs for " << iterations << " iterations\n";
  std::cout << "  Average: " << avg_us << " µs per molecule\n";
  if (tetrahedral + double_bonds > 0) {
    double us_per_center = avg_us / static_cast<double>(tetrahedral + double_bonds);
    std::cout << "  Per stereocenter: " << us_per_center << " µs\n";
  }
  std::cout << "\n";
}

int main(int argc, char** argv) {
  std::cout << "NewCIPLabeler Quick Benchmark\n";
  std::cout << "==============================\n\n";

  // If PDB files provided on command line, benchmark those
  if (argc > 1) {
    std::cout << "Benchmarking PDB files:\n";
    std::cout << "-----------------------\n\n";
    for (int i = 1; i < argc; ++i) {
      benchmarkPDB(argv[i], 100);  // 100 iterations for PDB files
    }
  }

  // Standard SMILES benchmarks
  std::cout << "Benchmarking SMILES:\n";
  std::cout << "--------------------\n\n";

  // Simple cases (should be fast - constitutional ranking)
  benchmark("Simple tetrahedral", "Br[C@H](Cl)F", 10000);
  benchmark("Simple phosphine", "C[P@](C1CCCC1)C1=CC=CC=C1", 5000);

  // Cases requiring shell expansion
  benchmark("Ethyl vs methyl", "CC[C@H](C)F", 5000);
  benchmark("Ring system", "C[C@H]1CCCCC1", 5000);

  // Complex symmetric cases (worst case - deep expansion)
  benchmark("Large symmetric ring", "C[C@H]1CCCCCCCCCC1", 1000);
  benchmark("Branched symmetric", "CC(C)[C@H](C(C)C)Br", 1000);

  // Multiple stereocenters
  benchmark("Two centers", "C[C@H](Cl)CC[C@H](Cl)C", 5000);
  benchmark("Three centers", "C[C@H](F)C[C@H](Cl)C[C@H](Br)C", 2000);

  return 0;
}
