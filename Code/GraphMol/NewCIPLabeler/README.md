# NewCIPLabeler - High-Performance CIP Implementation

## Overview

This is a high-performance implementation of CIP (Cahn-Ingold-Prelog) stereochemical labeling that uses **lazy shell expansion** to dramatically improve speed over the standard CIPLabeler.

**Current Status:** Phase 1 complete, Phase 2 partial (78% passing)

## Quick Start

```cpp
#include <GraphMol/NewCIPLabeler/NewCIPLabeler.h>

auto mol = "Br[C@H](Cl)F"_smiles;
RDKit::NewCIPLabeler::assignCIPLabels(*mol);

std::string code;
mol->getAtomWithIdx(1)->getPropIfPresent(RDKit::common_properties::_CIPCode, code);
// code == "S"
```

## Architecture

### Core Design Principles

1. **Lazy Shell Expansion** - Only expand molecular shells when needed to break ties
2. **Constitutional Fast Path** - ~90% of centers resolve with just atomic number/mass
3. **Safety Limits** - Hard limits prevent infinite loops during development
4. **Function-Based Design** - Prefer functions over classes, following RDKit patterns

### Key Optimization Strategy

Most real molecules resolve quickly:
- **Constitutional rules** (atomic number, isotope) resolve most cases instantly
- **First few shells** (1-3 bonds away) resolve nearly all remaining cases
- **Full graph expansion** rarely needed

### File Structure

```
Code/GraphMol/NewCIPLabeler/
├── NewCIPLabeler.h/cpp      # Public API (matches old CIPLabeler)
├── Descriptor.h/cpp         # Enum for R/S/E/Z/M/P + to_string()
├── Substituent.h/cpp        # Lightweight substituent + lazy shell expansion
├── Rules.h/cpp              # CIP priority rules (1a, 1b, 2)
├── Ranker.h/cpp             # Core ranking/comparison engine
├── Tetrahedral.h/cpp        # R/S tetrahedral labeling
├── DoubleBond.h             # E/Z double bond labeling (stub)
├── catch_tests.cpp          # Test suite
├── CMakeLists.txt           # Build configuration
├── PHASE2_STATUS.md         # Detailed Phase 2 status
└── README.md                # This file
```

## Test Results

### Phase 1: Foundation (Complete) ✅
**Status:** 33/33 assertions passing (100%)

**Coverage:**
- Constitutional ranking (atomic number)
- Isotope comparison
- Simple tetrahedral centers (Br, Cl, F, H)
- Descriptor enum
- Exception handling

**Example:**
```cpp
auto mol = "Br[C@H](Cl)F"_smiles;
assignCIPLabels(*mol);
// → "S" (correct)
```

### Phase 2: Shell Expansion (Partial) ⚠️
**Status:** 51/65 total assertions passing (78%)

**What Works:**
- Multi-shell expansion without hanging
- Safety limits (20 shells, 1000 atoms/shell)
- Lexicographic shell comparison
- Most tetrahedral centers

**What Doesn't:**
- clearProp() prevents re-labeling (mystery bug)
- Phosphine/arsine parity inverted
- Some molecules exceed shell limits

### Phase 3: Double Bonds (Not Started)
**Status:** Not implemented

**Planned:**
- E/Z double bond labeling
- Reuse shell expansion from Phase 2

## Known Issues

### 1. clearProp Mystery Bug (Critical) 🔴

**Problem:** Calling `atom->clearProp(_CIPCode)` causes ALL labeling to fail.

**Impact:** Cannot re-label atoms. SmilesToMol sets `_CIPCode` automatically, so fresh molecules work but re-labeling doesn't.

**Reproduction:**
```cpp
auto mol = "Br[C@H](Cl)F"_smiles;
auto atom = mol->getAtomWithIdx(1);
atom->clearProp(common_properties::_CIPCode);  // Clears existing label
assignCIPLabels(*mol);
// BUG: No label assigned! Expected "S"
```

**Status:** Root cause unknown. Workaround: Don't use clearProp().

### 2. Phosphine Parity Inversion (High Priority) 🟡

**Problem:** P/As centers get systematically inverted R/S labels.

**Impact:** All phosphine tests fail with opposite labels.

**Example:**
```cpp
auto mol = "C[P@](C1CCCC1)C1=CC=CC=C1"_smiles;
assignCIPLabels(*mol);
// Gets: "S"
// Expected: "R"
```

**Status:** Parity calculation or spatial order building has a bug. Regular carbon centers work correctly.

### 3. Shell Limit Too Low (Medium Priority) 🟡

**Problem:** `HARD_MAX_SHELLS = 20` is too conservative for some valid molecules.

**Impact:** Some molecules throw `MaxIterationsExceeded`.

**Workaround:** Increase limit to 50-100 (but increases risk on pathological cases).

## Safety Features

### Development Limits
```cpp
// Ranker.cpp
constexpr uint32_t HARD_MAX_SHELLS = 20;  // Max shell depth

// Substituent.cpp
constexpr size_t MAX_SHELL_SIZE = 1000;   // Max atoms per shell
if (parent_mult == 0) continue;            // Don't expand ring closures
```

### Test Timeout
```bash
TIMEOUT=30 ./run_newcip_tests.sh "[phase2]"
```

## API

### Main Entry Points

```cpp
namespace RDKit::NewCIPLabeler {

// Label all chiral centers
void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations = 0);

// Label specific atoms/bonds
void assignCIPLabels(ROMol &mol,
                     const boost::dynamic_bitset<> &atoms,
                     const boost::dynamic_bitset<> &bonds,
                     unsigned int maxRecursiveIterations = 0);

// Exception for convergence failures
class MaxIterationsExceeded : public std::runtime_error {};

}
```

### Output

Results stored in molecule properties:
- **Atoms:** `common_properties::_CIPCode` = "R", "S", "r", "s"
- **Bonds:** `common_properties::_CIPCode` = "E", "Z"
- **Molecule:** `common_properties::_CIPComputed` = true

## Building

```bash
# Build
/Users/dbn/builds/26-2/source/mmshare/build_tools/buildinger.sh --name 26-2 rdkit

# Run tests
./run_newcip_tests.sh "[phase1]"   # Phase 1 tests
./run_newcip_tests.sh "[phase2]"   # Phase 2 tests
./run_newcip_tests.sh "[newCIP]"   # All tests
```

## Performance

### Current Performance
- **Constitutional ranking:** Instant (no shell expansion)
- **1-3 shells:** Fast (typical case)
- **No benchmarking done yet** (planned for Phase 5)

### Expected Performance (vs old CIPLabeler)
- **Target:** 10x+ speedup on complex molecules
- **Mechanism:** Avoid full graph construction via lazy expansion
- **Trade-off:** Safety limits may reject some valid molecules

## Algorithm Details

### Ranking Flow

1. **Try Constitutional Ranking** (fast path)
   - Compare atomic numbers (Rule 1a)
   - Compare isotope masses (Rule 1b)
   - If all unique → done

2. **Incremental Shell Expansion**
   - Expand all substituents to shell 0 (root atom)
   - Apply CIP rules at depth 0
   - If not resolved → expand to shell 1
   - Repeat until resolved or max depth

3. **Shell Comparison**
   - Build descriptor: `(Z * 1000000 + mass * 1000 + multiplicity)`
   - Sort descriptors per shell (descending)
   - Lexicographic comparison across all shells
   - Higher descriptor = higher priority

4. **Parity Calculation** (tetrahedral)
   - Build CIP priority order (high → low)
   - Build spatial order (from bond iteration)
   - Count inversions → odd/even parity
   - Adjust chiral tag if odd parity
   - Assign: CCW → S, CW → R

### Shell Expansion Logic

```cpp
void Substituent::expandNextShell(const ROMol& mol, boost::dynamic_bitset<>& visited) {
  // Shell 0: root atom
  // Shell N: neighbors of shell N-1

  for (neighbor in previous_shell) {
    if (parent_mult == 0) continue;  // Skip duplicates (ring closures)

    for (bond in atom_bonds(neighbor)) {
      if (visited[nbr]) {
        add_as_duplicate(multiplicity = 0);  // Ring closure
      } else {
        add_new_atom(multiplicity = parent_mult * bond_order);
        visited.set(nbr);
      }
    }
  }
}
```

## Next Steps

### To Fix Phase 2 Issues

1. **clearProp Bug**
   - Trace through execution with debug logging
   - Check if SmilesToMol sets hidden state beyond `_CIPCode`
   - Consult RDKit developers about clearProp behavior

2. **Phosphine Parity**
   - Compare parity calculation for C vs P centers
   - Check if spatial order differs for P
   - Verify CIP priority order building

3. **Shell Limits**
   - Benchmark typical molecule complexity
   - Make `HARD_MAX_SHELLS` configurable
   - Add progressive timeout warnings

### To Complete Implementation

1. **Phase 3:** E/Z double bond labeling
2. **Phase 4:** Pseudo-asymmetry, atropisomers (M/P)
3. **Phase 5:** Performance benchmarks, optimization

### Testing Strategy

1. Import all tests from `Code/GraphMol/CIPLabeler/catch_tests.cpp` (1536 lines)
2. Run side-by-side comparison with old CIPLabeler
3. Benchmark performance on representative molecules

## References

**Algorithmic Paper:**
Hanson, R. M., et al. "Algorithmic Analysis of Cahn-Ingold-Prelog Rules of Stereochemistry."
*J. Chem. Inf. Model.* 2018, 58, 1755-1765.

**RDKit CIPLabeler (old implementation):**
`Code/GraphMol/CIPLabeler/`

## Contributors

- Implementation: Claude Sonnet 4.5
- Architecture: wonderboy + Claude
- Based on: RDKit CIPLabeler (Greg Landrum et al.)

## License

BSD License (matches RDKit)
