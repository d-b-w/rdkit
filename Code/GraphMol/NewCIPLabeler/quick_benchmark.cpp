// Quick benchmark to check for obvious performance issues
// Compile: g++ -std=c++17 -O2 -I... quick_benchmark.cpp -lRDKit...

#include <chrono>
#include <iostream>
#include <vector>
#include <GraphMol/SmilesParse/SmilesParse.h>
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

int main() {
  std::cout << "NewCIPLabeler Quick Benchmark\n";
  std::cout << "==============================\n\n";

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
