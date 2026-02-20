# High-Performance CIP Labeler Implementation Plan

## Context

The existing CIPLabeler in `Code/GraphMol/CIPLabeler/` is accurate but too slow for production use. Performance issues are fundamental to its design (digraph-based approach with full graph expansion). We need a new implementation that:

- **Matches accuracy** of existing implementation (pass all existing tests)
- **Dramatically improves performance** by avoiding full graph construction
- **Follows the algorithmic paper** rather than copying the old implementation
- **Maintains clean, modern C++ code** following RDKit patterns

The new implementation will live in `Code/GraphMol/NewCIPLabeler/` alongside the old one, allowing gradual migration and direct performance comparison.

## Architecture Overview

### Core Design Principles

1. **Lazy shell expansion**: Only expand molecular shells when needed, not upfront
2. **Constitutional fast path**: ~90% of centers resolve with just atomic number/mass (Rule 1)
3. **Batch processing**: Resolve simple centers first, complex ones last
4. **Function-based design**: Prefer functions over classes, static/unnamed namespaces for helpers
5. **Modern C++**: Use constexpr, lambdas, std::ranges, bitsets for performance

### Key Optimization Strategy

Most real molecules have CIP labels that resolve quickly:
- **Constitutional rules** (atomic number, isotope) resolve most cases
- **First few shells** (1-3 bonds away) resolve nearly all remaining cases
- **Full graph expansion** rarely needed

Implementation exploits this by:
- Attempting constitutional comparison first (no graph building)
- Expanding shells incrementally only when ties occur
- Processing easy centers before hard ones (allows label reuse)

## File Structure

```
Code/GraphMol/NewCIPLabeler/
├── NewCIPLabeler.h              # Public API (matches old CIPLabeler.h interface)
├── NewCIPLabeler.cpp            # Entry points, center collection, batching
├── Descriptor.h                 # Enum for R/S/E/Z/M/P/etc. + to_string()
├── Substituent.h                # Lightweight substituent representation
├── Substituent.cpp              # Shell expansion logic
├── Rules.h                      # CIP priority rules (1a, 1b, 2, etc.)
├── Rules.cpp                    # Rule implementations
├── Ranker.h                     # Core ranking/comparison engine
├── Ranker.cpp                   # Substituent comparison and sorting
├── Tetrahedral.h                # Tetrahedral (R/S) labeling
├── Tetrahedral.cpp
├── DoubleBond.h                 # E/Z double bond labeling
├── DoubleBond.cpp
├── catch_tests.cpp              # Test suite (duplicated from old CIPLabeler)
└── CMakeLists.txt               # Build configuration
```

## Public API (NewCIPLabeler.h)

Match existing interface exactly:

```cpp
namespace RDKit {
namespace NewCIPLabeler {

// Label all chiral centers in molecule
RDKIT_NEWCIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol,
    unsigned int maxRecursiveIterations = 0);

// Label specific atoms/bonds
RDKIT_NEWCIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol,
    const boost::dynamic_bitset<> &atoms,
    const boost::dynamic_bitset<> &bonds,
    unsigned int maxRecursiveIterations = 0);

// Exception for convergence failures
class RDKIT_NEWCIPLABELER_EXPORT MaxIterationsExceeded
    : public std::runtime_error {
  public:
    explicit MaxIterationsExceeded();
};

}  // namespace NewCIPLabeler
}  // namespace RDKit
```

Results stored in molecule properties:
- Atoms/bonds: `common_properties::_CIPCode` = "R"/"S"/"E"/"Z"/etc.
- Molecule: `common_properties::_CIPComputed` = true

## Data Structures

### Descriptor Enum (Descriptor.h)

```cpp
enum class Descriptor {
  NONE,      // Unspecified
  UNKNOWN,   // Cannot determine
  ns,        // Unspecified other

  // Tetrahedral
  R, S,      // Chiral
  r, s,      // Pseudo-chiral

  // Double bond
  E, Z,      // Entgegen/Zusammen
  seqCis, seqTrans,  // Sequential

  // Atropisomer
  M, P,      // Axial chirality
  m, p,      // Pseudo-axial

  // Coordination (future)
  SP_4, TBPY_5, OC_6
};

std::string to_string(const Descriptor &desc);
```

### Substituent Representation (Substituent.h)

Lightweight structure for lazy expansion:

```cpp
// Single shell of atoms at distance N from root
struct AtomShell {
  uint32_t distance;
  std::vector<const Atom*> atoms;
  std::vector<uint32_t> multiplicities;  // Bond order duplicates
};

// One substituent extending from a stereocenter
struct Substituent {
  const Atom* root_atom;           // Immediate neighbor
  const Bond* connecting_bond;     // Bond to center
  std::vector<AtomShell> shells;   // Expanded shells (lazy)

  // Expand one more shell
  void expandNextShell(const ROMol& mol, boost::dynamic_bitset<>& visited);
};

// Ranking result for a center
struct CenterRanking {
  bool is_unique;             // All substituents ranked?
  bool is_pseudo;             // Use r/s instead of R/S?
  std::vector<size_t> order;  // Sorted indices (low to high priority)
};
```

## Algorithm Flow

### 1. Top-Level Flow (NewCIPLabeler.cpp)

```cpp
void assignCIPLabels(ROMol& mol, unsigned int maxIters) {
  boost::dynamic_bitset<> all_atoms(mol.getNumAtoms());
  boost::dynamic_bitset<> all_bonds(mol.getNumBonds());
  all_atoms.set();
  all_bonds.set();
  assignCIPLabels(mol, all_atoms, all_bonds, maxIters);
}

void assignCIPLabels(ROMol& mol,
                     const boost::dynamic_bitset<>& atom_filter,
                     const boost::dynamic_bitset<>& bond_filter,
                     unsigned int maxIters) {
  // 1. Collect potential stereocenters
  auto atoms = findStereoAtoms(mol, atom_filter);
  auto bonds = findStereoBonds(mol, bond_filter);

  // 2. Sort by estimated complexity (simple first)
  auto centers = sortByComplexity(atoms, bonds, mol);

  // 3. Process each center
  for (const auto& center : centers) {
    labelCenter(mol, center, maxIters);
  }

  // 4. Mark as computed
  mol.setProp(common_properties::_CIPComputed, true);
}
```

### 2. Center Complexity Estimation

```cpp
namespace {

uint32_t estimateComplexity(const ROMol& mol, const Atom* atom) {
  uint32_t complexity = 0;

  // Ring membership increases complexity
  complexity += mol.getRingInfo()->numAtomRings(atom->getIdx()) * 10;

  // Heteroatoms in neighborhood increase complexity
  for (const auto& nbr : mol.atomNeighbors(atom)) {
    if (nbr->getAtomicNum() > 6) complexity += 5;
  }

  // Symmetry hints (same atomic number neighbors)
  std::unordered_map<int, int> z_counts;
  for (const auto& nbr : mol.atomNeighbors(atom)) {
    z_counts[nbr->getAtomicNum()]++;
  }
  for (const auto& [z, count] : z_counts) {
    if (count > 1) complexity += count * 3;
  }

  return complexity;
}

}  // namespace
```

### 3. Core Ranking Algorithm (Ranker.cpp)

```cpp
CenterRanking rankSubstituents(const ROMol& mol,
                               const Atom* center,
                               std::vector<Substituent>& subs,
                               uint32_t max_shells) {
  // Fast path: try constitutional rules only (no expansion)
  if (tryConstitutionalRanking(mol, subs)) {
    return buildRankingResult(subs, /*is_pseudo=*/false);
  }

  // Need deeper analysis - expand shells iteratively
  boost::dynamic_bitset<> visited(mol.getNumAtoms());
  visited.set(center->getIdx());  // Don't backtrack to center

  for (uint32_t shell = 0; shell < max_shells; ++shell) {
    // Expand all substituents to this shell
    for (auto& sub : subs) {
      if (sub.shells.size() <= shell) {
        sub.expandNextShell(mol, visited);
      }
    }

    // Apply CIP rules at this depth
    bool resolved = applyRankingRules(mol, subs, shell);
    if (resolved) {
      bool is_pseudo = checkPseudoAsymmetry(subs);
      return buildRankingResult(subs, is_pseudo);
    }
  }

  throw MaxIterationsExceeded();
}
```

### 4. Constitutional Fast Path (Rules.cpp)

```cpp
namespace {

// Rule 1a: Higher atomic number has priority
// Rule 1b: Higher isotope mass has priority
bool tryConstitutionalRanking(const ROMol& mol,
                              std::vector<Substituent>& subs) {
  // Build sorting keys: (atomic_num << 32) | isotope
  std::vector<std::pair<uint64_t, size_t>> keys;

  for (size_t i = 0; i < subs.size(); ++i) {
    const Atom* atom = subs[i].root_atom;

    // Handle implicit H (nullptr sentinel)
    if (atom == nullptr) {
      keys.emplace_back(1ULL << 32, i);  // H = Z=1, mass=1
      continue;
    }

    uint64_t z = atom->getAtomicNum();
    uint64_t mass = atom->getIsotope() ? atom->getIsotope() :
                    PeriodicTable::getTable()->getMostCommonIsotopeMass(z);

    uint64_t key = (z << 32) | mass;
    keys.emplace_back(key, i);
  }

  // Sort by key (higher = higher priority)
  std::sort(keys.begin(), keys.end(),
    [](const auto& a, const auto& b) { return a.first > b.first; });

  // Check for ties
  for (size_t i = 1; i < keys.size(); ++i) {
    if (keys[i].first == keys[i-1].first) {
      return false;  // Tie detected, need deeper analysis
    }
  }

  // All unique - assign ranks
  for (size_t rank = 0; rank < keys.size(); ++rank) {
    subs[keys[rank].second].final_rank = rank;
  }

  return true;
}

}  // namespace
```

### 5. Shell Expansion (Substituent.cpp)

```cpp
void Substituent::expandNextShell(const ROMol& mol,
                                  boost::dynamic_bitset<>& visited) {
  uint32_t new_shell_idx = static_cast<uint32_t>(shells.size());

  AtomShell new_shell{new_shell_idx, {}, {}};

  if (new_shell_idx == 0) {
    // Shell 0: just the root atom
    new_shell.atoms.push_back(root_atom);
    new_shell.multiplicities.push_back(1);
    visited.set(root_atom->getIdx());
  } else {
    // Expand from previous shell
    const auto& prev = shells[new_shell_idx - 1];

    for (size_t i = 0; i < prev.atoms.size(); ++i) {
      const Atom* atom = prev.atoms[i];
      uint32_t parent_mult = prev.multiplicities[i];

      for (const auto& bond : mol.atomBonds(atom)) {
        const Atom* nbr = bond->getOtherAtom(atom);

        if (visited[nbr->getIdx()]) {
          // Ring closure - add as duplicate node
          new_shell.atoms.push_back(nbr);
          new_shell.multiplicities.push_back(0);  // Zero indicates duplicate
          continue;
        }

        // New atom - add with bond order multiplicity
        int bond_order = static_cast<int>(bond->getBondTypeAsDouble());
        if (bond->getBondType() == Bond::DATIVE) {
          bond_order = 0;  // Dative bonds = 0
        }

        new_shell.atoms.push_back(nbr);
        new_shell.multiplicities.push_back(parent_mult * bond_order);
        visited.set(nbr->getIdx());
      }
    }
  }

  shells.push_back(std::move(new_shell));
}
```

### 6. Tetrahedral Labeling (Tetrahedral.cpp)

```cpp
namespace {

Descriptor computeTetrahedralDescriptor(const ROMol& mol,
                                        const Atom* center,
                                        const CenterRanking& ranking) {
  PRECONDITION(ranking.is_unique, "Cannot label without unique ranking");
  PRECONDITION(ranking.order.size() == 4, "Tetrahedral must have 4 substituents");

  // Get spatial arrangement from chiral tag
  auto tag = center->getChiralTag();
  if (tag != Atom::CHI_TETRAHEDRAL_CW && tag != Atom::CHI_TETRAHEDRAL_CCW) {
    return Descriptor::NONE;
  }

  // Compute parity based on ranking permutation
  // (Details depend on how RDKit stores stereo - see old tests)
  int parity = computeParity(mol, center, ranking.order);

  if (parity == 1) {
    return ranking.is_pseudo ? Descriptor::r : Descriptor::R;
  } else if (parity == 2) {
    return ranking.is_pseudo ? Descriptor::s : Descriptor::S;
  }

  return Descriptor::NONE;
}

}  // namespace

void labelTetrahedralCenter(ROMol& mol, Atom* center, uint32_t max_iters) {
  // Collect 4 neighbors (including implicit H if present)
  std::vector<Substituent> subs;

  for (const auto& bond : mol.atomBonds(center)) {
    const Atom* nbr = bond->getOtherAtom(center);
    subs.emplace_back(nbr, bond, std::vector<AtomShell>{});
  }

  // Add implicit H if needed
  if (center->getTotalNumHs() > 0 && subs.size() == 3) {
    subs.emplace_back(nullptr, nullptr, std::vector<AtomShell>{});
  }

  if (subs.size() != 4) {
    return;  // Not tetrahedral
  }

  // Rank substituents
  CenterRanking ranking = rankSubstituents(mol, center, subs, max_iters);

  if (!ranking.is_unique) {
    return;  // Cannot determine label
  }

  // Compute and assign descriptor
  Descriptor desc = computeTetrahedralDescriptor(mol, center, ranking);
  if (desc != Descriptor::NONE) {
    center->setProp(common_properties::_CIPCode, to_string(desc));
  }
}
```

## Implementation Phases (TDD)

### Phase 1: Foundation (1-2 days)
**Goal**: Minimal working implementation for simple tetrahedral cases

**Files**:
- `Descriptor.h` - enum definition
- `Substituent.h/.cpp` - data structures, shell expansion
- `Rules.h/.cpp` - constitutional comparison only
- `Ranker.h/.cpp` - ranking with constitutional fast path
- `Tetrahedral.h/.cpp` - basic tetrahedral labeling
- `NewCIPLabeler.h/.cpp` - public API stubs
- `catch_tests.cpp` - first 5 tests

**Tests**:
```cpp
TEST_CASE("Constitutional - atomic number", "[newCIP]") {
  auto mol = "Br[C@H](Cl)F"_smiles;  // Br > Cl > F > H
  NewCIPLabeler::assignCIPLabels(*mol);
  CHECK(mol->getAtomWithIdx(1)->getProp<std::string>(_CIPCode) == "S");
}

TEST_CASE("Constitutional - isotopes", "[newCIP]") {
  auto mol = "[2H][C@H](C)F"_smiles;  // F > C > D > H
  NewCIPLabeler::assignCIPLabels(*mol);
  // Check label...
}
```

**Acceptance**: 5+ simple tetrahedral tests passing

### Phase 2: Shell Expansion (2-3 days)
**Goal**: Handle cases requiring 1-3 shell expansions

**Enhancements**:
- Complete shell expansion in `Substituent.cpp`
- Full CIP rule implementation in `Rules.cpp`
- Lexicographic shell comparison
- Ring duplicate node handling

**Tests**:
```cpp
TEST_CASE("One shell expansion", "[newCIP]") {
  auto mol = "CC[C@H](C)CC"_smiles;  // Need to look at second carbon
  // ...
}

TEST_CASE("Ring systems", "[newCIP]") {
  auto mol = "C1CC[C@H](Br)CC1"_smiles;
  // ...
}
```

**Acceptance**: 20+ tetrahedral tests passing, including multi-shell cases

### Phase 3: Double Bonds (1-2 days)
**Goal**: E/Z labeling for double bonds

**Files**:
- `DoubleBond.h/.cpp` - double bond labeling

**Tests**:
```cpp
TEST_CASE("Simple E/Z", "[newCIP]") {
  auto mol = "F/C=C/Cl"_smiles;
  NewCIPLabeler::assignCIPLabels(*mol);
  CHECK(mol->getBondWithIdx(1)->getProp<std::string>(_CIPCode) == "E");
}
```

**Acceptance**: Double bond tests from old test suite passing

### Phase 4: Advanced Features (2-3 days)
**Goal**: Pseudo-asymmetry, atropisomers, edge cases

**Files**:
- `Atropisomer.h/.cpp` - M/P labeling
- Enhanced pseudo-asymmetry detection in `Ranker.cpp`

**Tests**: Import all tests from `Code/GraphMol/CIPLabeler/catch_tests.cpp`

**Acceptance**: All 1536 lines of tests passing

### Phase 5: Performance Optimization (1-2 days)
**Goal**: Ensure performance is dramatically better than old implementation

**Optimizations**:
- Batch processing by complexity
- Memoization of comparisons
- Benchmark tests

**Tests**:
```cpp
TEST_CASE("Performance - complex symmetric molecule", "[.benchmark]") {
  // Large molecule from old test suite
  auto start = std::chrono::high_resolution_clock::now();
  NewCIPLabeler::assignCIPLabels(*mol, 1000000);
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  // Compare to old implementation
}
```

**Acceptance**: 10x+ speedup on complex molecules

## Critical Files

### For Interface Reference
- `/Users/dbn/builds/26-2/source/rdkit/Code/GraphMol/CIPLabeler/CIPLabeler.h` - Public API to match

### For Test Duplication
- `/Users/dbn/builds/26-2/source/rdkit/Code/GraphMol/CIPLabeler/catch_tests.cpp` - 1536 lines of comprehensive tests

### For Style Reference
- `/Users/dbn/builds/26-2/source/rdkit/Code/GraphMol/Canon.cpp` - Function-based design, unnamed namespaces
- `/Users/dbn/builds/26-2/source/rdkit/Code/GraphMol/MolOps/FindRings.cpp` - Bitset usage, modern C++
- `/Users/dbn/builds/26-2/source/rdkit/Code/GraphMol/RDKitBase.h` - ROMol/Atom/Bond interfaces

### For Algorithm Reference
- `/Users/dbn/Downloads/algorithmic-analysis-of-cahn-ingold-prelog-rules-of-stereochemistry-proposals-for-revised-rules-and-a-guide-for-machine-implementation.pdf` - Authoritative CIP rules

## Verification Strategy

### During Development
1. **Each commit**: Run tests for current phase
2. **Each phase**: Ensure no regressions in previous tests
3. **Use git amend**: If implementation fails, update commit message explaining what was wrong

### Final Verification
1. **Accuracy**: All tests from old CIPLabeler pass
2. **Performance**: Benchmark on representative molecules, compare to old implementation
3. **Integration**: Can build with premake.py/CMake
4. **Documentation**: Code comments use terminology from algorithmic paper

### Test Commands
```bash
# Build
cd /Users/dbn/builds/26-2/source/rdkit
mkdir build && cd build
cmake ..
make NewCIPLabelerTests

# Run tests
./Code/GraphMol/NewCIPLabeler/NewCIPLabelerTests

# Benchmark (optional)
./Code/GraphMol/NewCIPLabeler/NewCIPLabelerTests "[.benchmark]"
```

## Success Criteria

1. ✓ All existing CIPLabeler tests pass
2. ✓ Performance 10x+ better on complex molecules
3. ✓ Code follows RDKit style (functions, unnamed namespaces, modern C++)
4. ✓ Incremental git commits document development process
5. ✓ Can coexist with old CIPLabeler for comparison
