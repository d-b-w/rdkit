# NewCIPLabeler Benchmarking Tools

## Quick Benchmark (quick_benchmark.cpp)

Runs standard SMILES-based benchmarks and optional PDB files.

### Build
```bash
cd /Users/dbn/builds/26-2/source/rdkit
./build_and_test.sh
```

This builds `benchmarkNewCIPLabeler` executable in the build directory.

### Usage

**Run standard benchmarks:**
```bash
cd /Users/dbn/builds/26-2/build/rdkit/build
./Code/GraphMol/NewCIPLabeler/benchmarkNewCIPLabeler
```

**Benchmark PDB files:**
```bash
./Code/GraphMol/NewCIPLabeler/benchmarkNewCIPLabeler /path/to/protein.pdb /path/to/ligand.pdb
```

### Output

For each test case, shows:
- Total time (µs)
- Average time per molecule (µs)
- Throughput (molecules/sec)
- For PDB: time per stereocenter (µs)

## PDB Test Tool (test_pdb.cpp)

Detailed analysis of a single PDB file.

### Build
Same as above - builds `testNewCIPLabeler_pdb` executable.

### Usage

```bash
cd /Users/dbn/builds/26-2/build/rdkit/build
./Code/GraphMol/NewCIPLabeler/testNewCIPLabeler_pdb /path/to/structure.pdb
```

### Output

Provides detailed information:
- Parse time
- Molecule statistics (atoms, bonds)
- Stereocenter counts (tetrahedral, double bonds)
- CIP labeling time
- Number of labeled centers
- First 10 labels as examples
- Performance metrics

## Example PDB Files to Test

### Small molecules with stereocenters:
- Drug-like molecules (100-500 atoms)
- Expected: 1-10 stereocenters
- Should be very fast (<100 µs)

### Medium proteins:
- Small peptides (1000-2000 atoms)
- Expected: 0-5 stereocenters (usually none)
- Tests scalability

### Large biomolecules:
- Proteins (5000+ atoms)
- Complex structures with many rings
- Tests worst-case performance

## Performance Expectations

Based on analysis in PERFORMANCE_ANALYSIS.md:

**Fast cases** (constitutional ranking):
- 5-50 µs per stereocenter
- ~90% of real molecules

**Medium cases** (1-3 shell expansion):
- 50-200 µs per stereocenter
- ~9% of real molecules

**Slow cases** (deep symmetric structures):
- 200-1000+ µs per stereocenter
- ~1% of real molecules
- Affected by O(N²) bottleneck (see PERFORMANCE_ANALYSIS.md)
