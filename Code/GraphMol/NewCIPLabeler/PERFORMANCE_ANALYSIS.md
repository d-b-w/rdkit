# NewCIPLabeler Performance Analysis

## Critical Bottlenecks Found

### 1. **Quadratic Key Rebuilding** (CRITICAL - Rules.cpp:113-164)

**Problem**: `applyRankingRules()` rebuilds comparison keys from scratch at each shell depth.

**Current behavior**:
```cpp
for (uint32_t shell = 0; shell < max_iter; ++shell) {
  applyRankingRules(mol, subs, shell);  // Rebuilds ALL shells 0..shell
}
```

At depth 0: compute shell 0
At depth 1: recompute shell 0 + compute shell 1
At depth 2: recompute shells 0, 1 + compute shell 2
...

**Complexity**: O(depth²) instead of O(depth)

**Impact**: For molecules requiring 10 shell expansions, we do 10x more work than needed

**Fix**: Store comparison keys in Substituent struct, append incrementally:
```cpp
struct Substituent {
  std::vector<int> comparison_key;  // Build incrementally
  // ...
};

// In applyRankingRules, only process NEW shell:
void appendShellToKey(Substituent& sub, uint32_t new_shell_depth) {
  const auto& shell = sub.shells[new_shell_depth];
  std::vector<int> shell_descriptors;
  // ... build descriptors for THIS shell only
  std::sort(shell_descriptors.rbegin(), shell_descriptors.rend());
  sub.comparison_key.insert(sub.comparison_key.end(),
                             shell_descriptors.begin(),
                             shell_descriptors.end());
}
```

**Estimated speedup**: 5-10x for molecules requiring deep expansion (rings, symmetric structures)

---

### 2. **Repeated Periodic Table Lookups** (HIGH - Rules.cpp:76-78, 147-149)

**Problem**: `getMostCommonIsotopeMass(z)` called repeatedly for same atomic numbers

**Current behavior**:
```cpp
// Called for every atom in every shell in every substituent
auto* pt = PeriodicTable::getTable();
unsigned int mass = atom->getIsotope() ? atom->getIsotope() :
    static_cast<unsigned int>(pt->getMostCommonIsotopeMass(z));
```

For a molecule with 100 carbon atoms across all shells, we look up carbon's mass 100 times.

**Impact**: Moderate - depends on PeriodicTable implementation (likely has internal cache)

**Fix**: Cache common masses in a static array:
```cpp
namespace {
constexpr int MAX_ATOMIC_NUM = 118;
std::array<unsigned int, MAX_ATOMIC_NUM + 1> common_masses;
std::once_flag masses_initialized;

void initCommonMasses() {
  auto* pt = PeriodicTable::getTable();
  for (int z = 1; z <= MAX_ATOMIC_NUM; ++z) {
    common_masses[z] = static_cast<unsigned int>(
        pt->getMostCommonIsotopeMass(z));
  }
}

unsigned int getCommonMass(int z) {
  std::call_once(masses_initialized, initCommonMasses);
  return common_masses[z];
}
}
```

**Estimated speedup**: 1.2-2x for molecules with many atoms

---

### 3. **Bond Order Switch Statement** (LOW - Substituent.cpp:76-97)

**Problem**: Manual switch statement instead of using RDKit's built-in method

**Current code**:
```cpp
int bond_order = 1;
switch (bond->getBondType()) {
  case Bond::SINGLE: bond_order = 1; break;
  case Bond::DOUBLE: bond_order = 2; break;
  case Bond::TRIPLE: bond_order = 3; break;
  // ... 7 cases total
}
```

**Better**:
```cpp
int bond_order = static_cast<int>(bond->getBondTypeAsDouble());
if (bond->getBondType() == Bond::DATIVE) {
  bond_order = 0;  // Special case for dative
}
```

**Impact**: Negligible - compiler likely optimizes switch to jump table

**Estimated speedup**: <1.05x

---

### 4. **Vector Allocations in Hot Loop** (MEDIUM - Rules.cpp:118, 130)

**Problem**: Allocating vectors inside nested loops

**Current code**:
```cpp
for (size_t i = 0; i < subs.size(); ++i) {
  std::vector<int> key;  // Allocation
  for (uint32_t depth = 0; depth <= shell_depth; ++depth) {
    std::vector<int> shell_descriptors;  // Allocation
    // ...
  }
}
```

**Fix**: Reuse buffers or reserve capacity upfront
```cpp
std::vector<int> key;
std::vector<int> shell_descriptors;
for (size_t i = 0; i < subs.size(); ++i) {
  key.clear();
  key.reserve(estimated_total_atoms);
  for (uint32_t depth = 0; depth <= shell_depth; ++depth) {
    shell_descriptors.clear();
    shell_descriptors.reserve(estimated_shell_size);
    // ...
  }
}
```

**Impact**: Moderate for molecules with many stereocenters

**Estimated speedup**: 1.2-1.5x

---

## Non-Critical Observations

### Good Practices Already Used ✓
- Bitset for visited tracking (O(1) lookups)
- Reserve() called on vectors where size is known
- Move semantics used (std::move)
- Early termination (constitutional fast path)
- Safety limits (MAX_SHELL_SIZE, HARD_MAX_SHELLS)

### Minor Improvements
- Lines 69-75 in buildRankingResult: std::map for rank counting
  - Could use std::unordered_map for O(1) instead of O(log n)
  - Impact: negligible (only 4 substituents per center)

---

## Recommended Priority

1. **Fix #1 (Incremental keys)**: Critical, easy fix, 5-10x speedup
2. **Fix #2 (Cache masses)**: High value, medium effort, 1.2-2x speedup
3. **Fix #4 (Buffer reuse)**: Medium value, easy fix, 1.2-1.5x speedup
4. **Fix #3 (Bond order)**: Low value, trivial fix, <5% speedup

**Total estimated speedup**: 8-30x for complex molecules

---

## Testing Strategy

Before optimizing:
1. Add benchmark test cases (large symmetric molecules, deep rings)
2. Profile with actual timing measurements
3. Compare NewCIPLabeler vs old CIPLabeler on benchmark set

After each optimization:
1. Run full test suite (ensure correctness)
2. Measure speedup on benchmark cases
3. Commit each optimization separately for bisectability
