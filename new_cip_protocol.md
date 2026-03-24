# CIP Rules (Proposed Version from Paper)

## Application Order

1. Build complete molecular digraph
2. Apply Rules 1a, 1b, 2 exhaustively to entire digraph
3. **Assign all auxiliary descriptors** (R/S/E/Z/seqCis/seqTrans/M/P/r/s/m/p/etc) to all nodes
4. Apply Rule 3: **Use** E/Z descriptors for ranking (Z > E)
5. Apply Rule 4a: **Use** descriptor types for ranking (chiral > pseudo-chiral > achiral)
6. Apply Rule 4b: **Use** like/unlike pairs for ranking
7. Apply Rule 5 if needed for pseudo-asymmetry detection

**Note on Auxiliary Descriptors**: These are the stereochemical labels assigned to nodes in the digraph: R, S (tetrahedral chirality), E, Z (double bond geometry), M, P (axial chirality), r, s (pseudo-tetrahedral), m, p (pseudo-axial), seqCis, seqTrans (sequential double bonds). They must be computed for ALL nodes before Rules 3-5 can use them for ranking.

## The Rules

**Rule 1a (proposed)**: Higher atomic number precedes lower atomic number.

**Rule 1b (proposed)**: Lower root distance precedes higher root distance.
- **Ring-closure duplicates**: distance = sphere of the duplicated atom
- **Multiple-bond duplicates**: distance = sphere of the atom the duplicate **attaches to** (NOT the duplicated atom)
- **All other atoms**: distance = sphere of the atom itself

**Important**: For multiple bonds, the duplicate node gets a **closer** distance than the real atom it duplicates!

Example: `Root -- C = C -- OH`
- First C is at distance 1 from root
- Second C is at distance 2 from root
- Double bond creates a duplicate of second C
- That duplicate gets distance **1** (same as first C, which it attaches to)
- So the duplicate (distance 1) ranks **higher** than the real atom (distance 2)

This avoids Kekulé structure ambiguities by treating the duplicate as being at the same sphere as the atom bearing the multiple bond.

**Rule 2 (proposed)**: Higher mass precedes lower mass.
- Duplicate nodes: mass = 0
- Atom with explicit isotope: exact isotopic mass
- All other atoms: atomic weight from periodic table
- Exception: For ¹⁶O, ⁵²Cr, ⁹⁶Mo, ¹⁷⁵Lu, subtract 0.1 from mass number when comparing

**Rule 3**: Z > E (nodes with Z or seqCis precede nodes with E or seqTrans)

**Rule 4a**: Chiral > pseudo-chiral > achiral.
- (R or S) > (r or s) > none
- (M or P) > (m or p) > none
- (seqCis or seqTrans) > (seqcis or seqtrans) > none

**Rule 4b**: Like descriptor pairs > unlike descriptor pairs.
- Normalize all descriptors to R or S
- Build priority key for each node: (path from root) + (Rule 4a priorities) + (like/unlike to reference)
- Ignore pseudo-asymmetric descriptors (r, s, m, p, seqcis, seqtrans) in comparison

**Rule 5**: Invert each descriptor and re-rank. If inversion reverses the order, center is pseudo-asymmetric (r/s instead of R/S).

## Implementation Details

### Multi-Phase Processing

Process ALL stereocenters through each phase before proceeding to the next:

**Phase 1: Build Complete Digraphs (Rules 1a, 1b, 2)**
- For each stereocenter, build its complete digraph by expanding shells incrementally
- During expansion, apply Rules 1a, 1b, 2 to rank nodes at each shell
- Incremental building is valid because these rules use only local information (atomic number, root distance, mass)
- Stop expanding when all substituents are uniquely ranked OR max depth reached
- If not uniquely ranked after Rules 1-2, the digraph remains for later phases

**Phase 2: Assign Auxiliary Descriptors (Outside-In)**

For each stereocenter being labeled:
1. Use the complete digraph from Phase 1
2. Assign auxiliary descriptors starting from the **outermost sphere** and working **inward** toward the root
3. When assigning descriptor to an auxiliary center at sphere N, all centers in spheres N+1, N+2, ... already have descriptors
4. Each auxiliary center only depends on spheres further from the root - never on anything between it and the root
5. The path from any auxiliary center back to root is unique (guaranteed by digraph) and can always be ranked by Rule 1a (atomic number)
6. Apply full CIP algorithm (Rules 1-5) recursively when computing each auxiliary descriptor

**Why this works:**
- No circular dependencies because we work from far-to-near
- Auxiliary center B's descriptor depends only on what's beyond B (away from root)
- The unique path back to root provides a tiebreaker (Rule 1a)
- Pseudo-asymmetry (r, s) is detected by Rule 5 during this recursive process

**Result:**
- All nodes in all digraphs have their auxiliary descriptors (R/S/E/Z/r/s/etc)
- These are **global properties** - once computed for an atom, that descriptor is used in all digraphs that include that atom
- All auxiliary descriptors must be assigned before Phase 3

**Phase 3: Apply Stereochemical Rules (Rules 3, 4a, 4b, 5)**
- For stereocenters not resolved in Phase 1, apply Rules 3-5 using the global auxiliary descriptors
- Rule 3: Use E/Z labels to rank nodes (Z > E)
- Rule 4a: Use descriptor types to rank (chiral > pseudo-chiral > achiral)
- Rule 4b: Use like/unlike pairs to rank
- Rule 5: Test for pseudo-asymmetry by inverting descriptors

### Key Points

1. **Incremental digraph building**: Valid for Rules 1a, 1b, 2 only (they use local information)
2. **Complete digraph required**: Rules 3-5 require the complete digraph with all auxiliary descriptors assigned
3. **Global auxiliary descriptors**: E/Z, R/S, M/P are molecular properties, not computed per-digraph
4. **Phase-by-phase, not atom-by-atom**: Process all atoms through Phase 1, then all through Phase 2, then all through Phase 3
5. **Duplicate nodes**: Ring closures and multiple bonds create duplicate nodes with mass=0
6. **Like/unlike**: Rule 4b compares descriptor sequences as binary strings (like=1, unlike=0)
