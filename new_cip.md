We'll plan C++ code to calculate Cahn−Ingold−Prelog (CIP) labels for molecules.

The implementation should follow the all rules (including the "proposed" rules) from "Algorithmic Analysis of Cahn−Ingold−Prelog Rules of Stereochemistry: Proposals for Revised Rules and a Guide for Machine Implementation" - see /Users/dbn/Downloads/algorithmic-analysis-of-cahn-ingold-prelog-rules-of-stereochemistry-proposals-for-revised-rules-and-a-guide-for-machine-implementation.pdf. Comments, variables, and functions should use terminology from the paper. We will be focused on performance as well as accuracy, so may take some short-cuts.

There is an implementation in Code/GraphMol/CIPLabeler. It's super slow, and the performance is irredemably tied to the implemetation. You are only allowed to read the interface (CIPLabeler.h and CIPLabeler.cpp) and test files (*catch*) from that directory! We should match the interface from Code/GraphMol/CIPLabeler/CIPLabeler.h - a few functions to label a molecule or subset thereof. The tests for CIPLabeler are great, we can duplicate them.

Internally, the implementation should sort/rank the substituents of each center. It should denote whether the sorted order is achiral, chiral, or pseudochiral. The implementation must be flexible enough to deal with both atom- and bond- centered stereochemistry.

CIP Chirality is a thorny problem to calculate efficiently and accurately. Our implementation should prioritize accuracy and maintainability, but performance is also very important. As we think of potential optimizations, you can take advantage of these assumptions.
* All potential stereo centers and bonds will be marked
* Most potential stereo centers will have a CIP label
* The CIP of most real molecules is resolved in the 3 few shells - we don't need to build the full graph to resolve
* Most CIP ranking are resolved with the "constitutional" rules. The full rule set is very expensive
* Most molecules have no isotopes - it's worth checking before even running rule 2.
* Proteins, cubane, and stereo-dependant stereo molecules are good performance tests.

The style of the code should follow the rest of Code/GraphMol/, with particular attention to Canon, MolOps, FindRings, ROMol. Most logic should be implemented in functions - don't use virtual inheritance. Functions/classes/struct should be local to the compilation unit when possible (static or unnamed namespace). Use constexpr when possible. Values should be const and locally scoped. Use C++20/17/14 features, including lamdbas, std::span, std::ranges. Stack/queue implementations are better than recursion. Typed functions are better than auto functions > template functions > virtual.

Development should follow a TDD approach, running successively more rigorous tests. We'll start with accuracy tests and move on to performance. All changes should be added as git commits. If an attempt fails, update the previous git commit message with an explanation of what is wrong with it.

## Safety
This code will include graph searching - it's easy to accidentally write infinite loops in graph search problems. It's OK to write hard depth limits on searches until the implementation is complete. It is also important to time-out test execution externally.

1. For which steps is valid to incrementally build the digraph?
2. For rule 3 - is E/Z a _global_ E/Z, or E/Z from the perspective of the incoming graph search?
3. For rule 4a - The rules before 4a don't talk about pseudo-chiral. How can there be a pseudo-chiral node at this point? (I beleive there is, I just don't understand why)
