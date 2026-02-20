# New CIP Implementation - Appendix

This document clarifies ambiguities in `new_cip.md` and provides additional context for implementation.

## Code Style Guidelines

Based on Canon, MolOps, FindRings, and ROMol:

### Code Organization
- **Functions over classes**: Most logic in free functions within namespaces (e.g., `RDKit::Canon::`, `RDKit::MolOps::`)
- **Unnamed namespace for file-local helpers**: `namespace { ... }` for functions only used in one file
- **`details` namespace for internal helpers**: `namespace details { ... }` for implementation details that need to be testable
- **Static when appropriate**: For truly private compilation unit functions

### C++ Features to Use
- **Lambdas extensively**: Both inline and for local algorithms
- **`auto` for type deduction**: Especially in range-for loops
- **Range-based for loops**: `for (auto atom : mol.atoms()) { ... }`
- **`constexpr` for constants**: `constexpr size_t MAX_BFSQ_SIZE = 200000;`
- **`std::ranges`**: Use where appropriate (C++20)
- **`std::span`**: Use for array views (C++20)
- **Structured bindings**: `auto [start, end] = ...;`

### Const Correctness
- **Extensive use of `const`**: Parameters, return values, member functions
- **Const references for inputs**: `const ROMol &mol`
- **Values are const and locally scoped**: Declare variables close to use

### Naming Conventions
- **CamelCase for types**: `ROMol`, `MolGraph`, `DigraphNode`
- **camelCase for functions/vars**: `getAtomWithIdx()`, `atomDegrees`, `computeRank()`
- **UPPER_CASE for typedefs/consts**: `RINGINVAR_SET`, `MAX_BFSQ_SIZE`
- **Member prefixes**: `d_` for data members, `dp_` for data pointers

### Documentation
- **Doxygen style**: `//!` and `/*! ... */`
- **Brief imperative form**: "Finds smallest rings" not "This function finds..."
- **Focus on "why" not "what"**: Comments explain rationale, not obvious code
- **No temporal references**: No "recently refactored" etc.

### Error Handling
- **PRECONDITION/CHECK_INVARIANT**: For assertions
- **Exceptions**: `ValueErrorException`, etc.
- **Fail early**: Return immediately on error conditions
- **No warnings in library code**: Either error or be silent

## Clarifications on Meta Plan Ambiguities

### 1. "Don't read it!" vs. "Match the interface"

**Resolution**: Read the existing interface (CIPLabeler.h) to understand the public API signature, but don't read the implementation (CIPLabeler.cpp). This allows matching the interface without being influenced by the slow implementation.

Steps:
1. Read `Code/GraphMol/CIPLabeler/CIPLabeler.h` for public API
2. Read test files to understand expected behavior
3. Do NOT read `Code/GraphMol/CIPLabeler/*.cpp` implementation files

### 2. "Don't use the paper as a recipe"

**Interpretation**: The paper provides the **rules** and **correctness criteria**, but not necessarily the optimal **algorithm**.

- **Follow**: The sequence rules (1a, 1b, 2, 3, 4a, 4b, 4c, 5, proposed 6)
- **Follow**: Terminology (digraph, duplicate node, auxiliary descriptor, etc.)
- **Follow**: Correctness requirements (what should be R vs S, etc.)
- **Adapt**: Implementation strategies for performance
- **Optimize**: Early termination, caching, incremental graph building

The paper is the **specification**, not the **recipe**.

### 3. Performance Optimization Strategy

The plan lists assumptions but no concrete strategy. Here's the approach:

**Leverage these assumptions:**
- All potential stereo centers and bonds will be marked (by preprocessing)
- Most potential stereo centers will have a CIP label (optimize common case)
- Most CIP labels resolve in first few shells (2-3 shells for typical molecules)
- Most rankings resolve with constitutional rules (Rules 1a, 1b, 2)

**Optimization strategies:**
1. **Incremental graph expansion**: Build digraph shell-by-shell, stop when decision made
2. **Lazy auxiliary descriptor calculation**: Only compute when needed by Rule 3+
3. **Early termination**: Stop at first rule that distinguishes ligands
4. **Caching**: Cache digraph nodes/rankings between centers when structure allows
5. **Shell depth limits**: Default max depth ~5-7 shells, configurable for edge cases

**Performance testing approach:**
- Start with correctness on validation suite
- Profile on realistic molecule sets
- Optimize hot paths identified by profiling

### 4. Who Marks Stereo Centers?

**Assumption**: The input molecule will have potential stereo centers marked via RDKit's existing `assignStereochemistry()` or similar preprocessing.

The new CIP labeler should:
- Accept a molecule with marked centers
- Compute correct CIP labels for marked centers
- Not be responsible for identifying which atoms/bonds are stereogenic

### 5. Data Structures

**Key data structures needed:**

```cpp
// Digraph node (real atom, duplicate for multiple bond, or ring closure)
struct DigraphNode {
    enum Type { REAL_ATOM, MULTIPLE_BOND_DUPLICATE, RING_CLOSURE_DUPLICATE };
    Type type;
    unsigned int atomIdx;  // index in molecule
    unsigned int sphere;   // distance from root
    double mass;           // for Rule 2 (0 for duplicates)
    // ... auxiliary descriptor if needed
};

// Digraph for a single center
struct Digraph {
    const ROMol* mol;
    unsigned int rootIdx;
    std::vector<DigraphNode> nodes;
    // adjacency or similar structure
};

// Ranking result for ligands
struct LigandRanking {
    std::vector<unsigned int> ligandOrder;  // indices into ligand list
    unsigned int decidingRule;  // which rule made the decision
    bool decided;
};
```

Prefer:
- `std::vector` for contiguous storage
- `std::span` for views into vectors
- `boost::dynamic_bitset` for sets of atoms/bonds (consistent with RDKit)

### 6. Git Workflow

**Revised workflow** (do NOT amend commits):

```bash
# Good commit messages during TDD
git commit -m "Add Rule 1a tests for simple cases"
git commit -m "Implement Rule 1a atomic number comparison"
git commit -m "Add test for mancude ring averaging"
git commit -m "Fix Rule 1a to handle mancude rings"
```

If an approach fails:
```bash
# Don't amend! Make a new commit explaining the issue
git commit -m "Revert incremental graph approach - insufficient for Rule 4b

The incremental expansion doesn't maintain enough state for
like/unlike comparisons in Rule 4b. Will need full digraph
generation before Rule 3 as paper specifies."
```

### 7. Test Sequencing

**Test progression** (aligned with TDD):

Phase 1 - Constitutional Rules:
1. Rule 1a: Simple atomic number cases (VS013-VS032)
2. Rule 1a: Mancude rings with averaging (VS032-VS033)
3. Rule 1b: Simple ring closures (VS171-VS174)
4. Rule 2: Isotopes (VS175-VS187)

Phase 2 - Stereochemistry-free:
5. Rule 3: Double bond configuration (VS188-VS195)

Phase 3 - Auxiliary Descriptors:
6. Rule 4a: Presence of stereogenicity (VS249-VS251)
7. Rule 4b: Like/unlike descriptors (VS196-VS263)
8. Rule 4c: Pseudoasymmetric priorities (VS273-VS279)
9. Rule 5: Enantiomorphic ligands (VS205-VS300)

Phase 4 - High Symmetry:
10. Rule 6: Spiro and axial cases (VS280-VS300)

Each phase builds on previous phases being correct.

### 8. Which Rule Versions to Implement?

**Implement the FULL PROPOSED RULE SET** from the paper, which includes 9 sequence rules:

**Rules with proposed modifications (use new versions):**

- **Rule 1b (proposed)**: Lower root distance precedes higher root distance, where root distance is defined:
  - (a) Ring-closure duplicates: sphere of the duplicated atom
  - (b) Multiple-bond duplicates: sphere of the atom to which duplicate is attached
  - (c) All other cases: sphere of the atom itself

- **Rule 2 (proposed)**: Higher mass precedes lower mass, where mass is defined as:
  - Duplicate nodes: 0
  - Atoms with isotope indicated: exact isotopic mass
  - All other cases: atomic weight

- **Rule 6 (proposed - NEW)**: An undifferentiated reference node has priority over any other undifferentiated node
  - Handles spiro compounds, axial symmetry cases
  - Apply after Rule 5 when 2+ pairs of identical ligands or 3-4 identical ligands remain

**Rules without modifications (use as specified in paper):**

- **Rule 1a**: Higher atomic number precedes lower (with mancude ring averaging)
- **Rule 3**: seqcis (Z) > seqtrans (E) > nonstereogenic
- **Rule 4a**: Chiral > pseudoasymmetric > nonstereogenic
- **Rule 4b**: Like descriptor pairs before unlike pairs
- **Rule 4c**: r > s, m > p
- **Rule 5**: R/M/seqCis > S/P/seqTrans (for enantiomorphic ligands)

**Rationale**: The paper demonstrates the proposed modifications are necessary for correctness on edge cases (mancude rings, isotopes, spiro compounds). The full rule set ensures complete coverage of the validation suite (VS001-VS300).

**Implementation note**: All 9 rules should be implemented. Do not implement BB-2013 versions alongside - use only the proposed rule set for simplicity and correctness.

### 9. Feature Prioritization

**MVP (Minimum Viable Product)** includes:
- **All 9 sequence rules**: 1a, 1b (proposed), 2 (proposed), 3, 4a, 4b, 4c, 5, 6 (proposed)
- **Tetrahedral centers**: Standard sp3 stereocenters
- **Double bonds (E/Z)**: Standard alkene stereochemistry
- **Spiro compounds**: Via Rule 6
- **Validation suite**: All VS001-VS300 passing where applicable

**Phase 2** additions (stereogenic unit types):
- Odd cumulenes (allenes, chirality axis)
- Atropisomers (restricted rotation)
- Even cumulenes (planar)

**Phase 3** (if needed - exotic cases):
- Helical chirality
- High-symmetry compounds beyond basic spiro
- In/out stereochemistry (bicyclooctanes, as noted in paper but not recommended)

**Prioritization strategy**:
- Implement all 9 rules for tetrahedral/double bond cases first (core MVP)
- Extend to other stereogenic unit types (allenes, atropisomers) only after core is working
- Don't over-engineer - the paper's validation suite defines the scope

### 10. Edge Cases and Error Handling

**Strategy**:
- Input validation: Check that stereochemistry is marked
- Digraph validation: Detect cycles that shouldn't exist
- Rule timeout: Maximum depth/iterations to prevent infinite loops
- Graceful degradation: Return "no descriptor" rather than crash on pathological cases

**Error cases**:
- Unmarked stereo center → skip (not our job to find them)
- Undecidable after Rule 5 → return achiral
- Graph too large/complex → configurable limits, warn and skip

### 11. Thread Safety

**Not a priority** for initial implementation. RDKit molecules are generally processed sequentially. Document that the implementation is not thread-safe without external synchronization.

If needed later: Make `const` member functions truly const (no internal mutation).

## Reading Guide for Existing CIPLabeler

**DO READ:**
- `Code/GraphMol/CIPLabeler/CIPLabeler.h` - Public interface
- `Code/GraphMol/CIPLabeler/test*.cpp` - Test cases and expected behavior
- Test molecules in validation suite

**DO NOT READ:**
- `Code/GraphMol/CIPLabeler/*.cpp` (implementation files)
- Any internal helper classes/functions

This preserves interface compatibility while avoiding implementation bias.

## References

- **Paper**: /Users/dbn/Downloads/algorithmic-analysis-of-cahn-ingold-prelog-rules-of-stereochemistry-proposals-for-revised-rules-and-a-guide-for-machine-implementation.pdf
- **RDKit Patterns**: Code/GraphMol/{Canon,MolOps,FindRings,ROMol}
- **Validation Suite**: https://cipvalidationsuite.github.io/ValidationSuite (VS001-VS300)
- **Blue Book 2013**: IUPAC Nomenclature of Organic Chemistry Chapter 9

## Glossary

- **Digraph**: Finite acyclic directed graph representing molecular structure for CIP analysis
- **Duplicate node**: Node added to digraph as copy of real atom (for multiple bonds or ring closures)
- **Duplicated atom**: Real atom that has been duplicated in the digraph
- **Auxiliary descriptor**: Temporary R/S/r/s assignment made during analysis of a specific center
- **Sphere**: Distance from root in digraph (shell number)
- **Mancude ring**: Ring having maximum number of noncumulative double bonds (e.g., benzene)
- **Ligand**: Substituent attached to stereogenic center
- **Constitutional**: Relating to connectivity/atomic properties (Rules 1-2)
- **Stereogenic**: Relating to spatial arrangement/stereochemistry (Rules 3-5)
