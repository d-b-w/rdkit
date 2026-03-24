/*
 *
 *
 *
 */

// not sure if this is needed
struct CIPMol {

    vector<Atom> atoms;
    vector<int> original_indices;

    CIPMol(ROMol& m) {
        // populate atoms?
        // reorder atoms by connectivity?
    }

};


void label(ROMol mol) {
    CIPMol cip_mol(mol);

    // find potential stereo
    for (auto a: mol.atoms()) {

    }
    //
    // do the consti
    //
    // sort

}

struct Atom {
int atomicNumber;
int isotope;
};

auto span(v, indices) {
    return std::span(v.begin() + indices[0], v.begin() + indices[1]);
}

void sort(v, indices, cmp)
{
    std::sort(v.begin() + indices[0], v.begin() + indices[1], cmp);
}


enum class CMP {
LESS, PSEUDO_LESS, EQUAL, PSEUDO_GREATER, GREATER
};
static constexpr auto EMPTY = std::limits<size_t>::max;

enum class BONDDESCRIPTOR {
    UNASSIGNED, NONE, Z, E
};
enum class ATOMDESCRIPTOR {
    UNASSIGNED, NONE, r, s, R, S
};


struct CIPAtom {
    CIPAtom(const CIPAtom&) = delete;
    CIPAtom(CIPAtom&&) = default;
    auto operator=(const CIPAtom&) = delete;
    auto operator=(CIPAtom&&) = default;

    uchar atomicNumber; // this may be fractional due to averaging
    double massNumber; // this may be fractional due to averaging
    size_t depth;
    BONDDESCRIPTOR bondDescriptor = BONDDESCRIPTOR::UNASSIGNED;
    ATOMDESCRIPTOR atomDescriptor = ATOMDESCRIPTOR::UNASSIGNED;
    // Descriptor descriptor = UNASSIGNED;
    // Descriptor auxDescriptor = UNASSIGNED;

    Atom* realAtom = nullptr;
    CIPAtom* sourceAtom = nullptr;
    CIPAtom* parent = nullptr;
    size_t idx = EMPTY;
    // char flags; // ring duplicate, bond duplicate, implicit hydrogen
    bool ringDuplicate = false;
    bool bondDuplicate = false;
    bool implicitHydrogen = false;
    bool expanded = false;

    // uchar shell_rank = 0;
    // uchar cross_rank = 0;

    std::span<CIPAtom> children;

    CIPAtom(CIPAtom* parent, CIPAtom& reference);
    CIPAtom(CIPAtom* parent, Atom* reference, bool bondDuplicate);

    bool isVisited(size_t i) const {
        // if it's not in a ring, this isn't possible
        auto a = this;
        while (a.idx != i && a.parent) {
            a = a.parent;
        }
        return a.idx == i;
    }

    unsigned int ringClosureDepth() const {
        // if it's not in a ring, this isn't possible
        if (realAtom) {
            auto ri = realAtom->getOwningMol()->getRingInfo();
            if (ri && ri->numAtomRings(idx) == 0) {
                return depth;
            }
        }

        auto a = this;
        const auto i = a->idx;
        while (a->idx != i && a->parent) {
            a = a->parent;
        }
        return a->depth;
    }

    size_t neighborCount() const {
        if (sourceAtom) {
            // based on an existing CIP DAG. use that count (+ the parent)
            return sourceAtom->children.size() + 1;
        } else {
            size_t count = 0u;
            count += realAtom->getTotalNumHs();
            for (auto b: realAtom->getOwningMol()->atomBonds()) {
                count += static_cast<size_t>(b->getBondTypeAsDouble());
            }
            return count;
        }
    }
    void expand(std::vector<CIPAtom>& atoms) {
        if (expanded) {
            return;
        }
        const auto initial_size = atoms.size();
        if (sourceAtom) {
            for (auto c: sourceAtom->children) {
                // duplicate it
                atoms.emplace_back({*this, c});
            }
        } else {
            for (auto b: realAtom->getOwningMol()->atomBonds()) {
                auto other = b->getOtherAtom(realAtom);
                auto reps = static_cast<size_t>(b->getBondTypeAsDouble());
                for (size_t i = 0; i < reps; ++i) {
                    atoms.emplace_back({*this, other, i != 0});
                }
            }
            for (size_t i = 0; i < realAtom->getTotalNumHs(); ++i) {
                atoms.emplace_back({*this, nullptr, i != 0});
            }
        }
        children = {atoms.begin() + initial_size, atoms.end()};
        expanded = true;
    }
};

CIPAtom::CIPAtom(CIPAtom* parent, CIPAtom& reference):
atomicNumber(reference.atomicNumber), massNumber(reference.massNumber), depth(parent->depth + 1),
realAtom(reference.realAtom), parent(parent), idx(reference.idx), ringDuplicate(reference.ringDuplicate),
bondDuplicate(reference.bondDuplicate), implicitHydrogen(reference.implicitHydrogen)
{

}
CIPAtom::CIPAtom(CIPAtom* parent, Atom* reference, bool bondDuplicate):
    depth(parent->depth + 1), realAtom(reference), parent(parent), bondDuplicate(bondDuplicate)
{
    if (reference) {
        atomicNumber = reference->getAtomicNum();
        massNumber = reference->getMass();
        if (bondDuplicate) {
            depth = parent->depth;
        } else {
            depth = ringClosureDepth();
        }
    } else {
        implicitHydrogen = true;
        atomicNumber = 1;
        massNumber = 1; // right?
    }
}

struct Shells;

struct ShellIterator
{
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = std::span<std::span<CIPAtom>>;
    using pointer = std::span<std::span<CIPAtom>>*;
    using reference = std::span<std::span<CIPAtom>>&;

    ShellIterator(Shells& shells, bool at_end)
        : shells(shells), at_end(at_end) {}

    ShellIterator& operator++() {
        if (at_end) {
            throw std::out_of_range("ShellIterator out of range");
        }
        ++current;
        if (shells.getShell(current) == nullptr) {
            at_end = true;
        }
        return *this;
    }

    ShellIterator operator++(int) {
        ShellIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator==(const ShellIterator& other) const {
        return current == other.current;
    }

    bool operator!=(const ShellIterator& other) const {
        return !(*this == other);
    }

    value_type& operator*() const {
        return *(shells.getShell(current));
    }

private:
    Shells& shells;
    size_t current = 0;
    bool at_end = false;
};

struct Shells
{
    Shells(const Atom* src, const Atom* ligand);

    // iterate over the shells out from something
    ShellIterator* begin() {return ShellIterator(this, false);}
    ShellIterator* end() {return ShellIterator(this, true);}
    std::span<std::span<CIPAtom>>* getShell(size_t idx);


    // allow re-use of the memory
    void assign(const Atom* src, const Atom* ligand);
    void clear()

private:
    // populate the next shell and fix spans
    bool reserveForNextShell();
    bool makeNextShell();

    std::vector<char> ever_visited;

    std::vector<CIPAtom> atoms; // all atoms
    std::vector<std::span<CIPAtom>> neighbor_groups; // groups of neighbors of some atom
    std::vector<std::span<std::span<CIPAtom>>> shells; // iterable of neighbor groups
};

Shells::Shells(const Atom* src, const Atom* ligand)
{
    assign(src, ligand);
}
void Shells::assign(const Atom* src, const Atom* ligand)
{
    clear();

    atoms.emplace_back(src, nullptr);
    atoms.emplace_back(ligand, src);
    ever_visited = {True, False};
    neighbor_groups.emplace_back(atoms.begin() + 1, 1);
    neighbor_groups_idx.emplace_back(1, 2);
    shells.emplace_back(neighbor_groups.begin(), 1);
    shells_idx.emplace_back(0, 1);
}
void Shells::clear()
{
    atoms.clear();
    ever_visited.clear()
    neighbor_groups.clear();
    neighbor_groups_idx.clear();
    shells.clear();
    shells_idx.clear();
}

std::span<std::span<CIPAtom>>* Shells::getShell(size_t idx)
{
    while (idx >= shells.size()) {
        if (!makeNextShell()) {
            return nullptr;
        }
    }
    return &shells[idx];
}

bool Shells::reserveForNextShell()
{
    // populate the next shell and fix spans
    size_t new_atom_count = 0;
    size_t newNeighbor_group_count = 0;
    for (auto& a : shells.back() | std::views::join) {
        const auto neighbors = a.neighborCount() - 1; // ignore parent
        newNeighbor_group_count += 1 ? neighbors != 0 : 1;
        new_atom_count += neighbors;
    }
    if (new_atom_count == 0) {
        return false;
    }
    const auto atoms_ptr = std::static_cast<void*>(atoms.data());
    atoms.reserve(atoms.size() + new_atom_count);
    const auto new_atoms_offset = std::static_cast<void*>(atoms.data()) - atoms_ptr;
    if (new_atoms_offset) {
        for (auto& a: atoms) {
            a.parent += new_atoms_offset;
        }
    }

    const auto neighbor_groups_ptr = std::static_cast<void*>(neighbor_groups.data());
    neighbor_groups.reserve(neighbor_groups.size() + newNeighbor_group_count);
    // reallocation, need to update span pointers
    if (new_atoms_offset != 0) {
        for (auto& g: neighbor_groups) {
            g.data_ += new_atoms_offset;
        }
    }
    const auto neighbor_groups_offset = std::static_cast<void*>(neighbor_groups.data());
    shells.reserve(shells.size() + newNeighbor_group_count);
    // reallocation, need to update span pointers
    if (neighbor_groups_offset != 0) {
        for (auto& s: shells) {
            s.data_ += neighbor_groups_offset;
        }
    }

    return true;
}

bool Shells::makeNextShell()
{
    // do all allocation ahead to allow use of spans for data
    if (!reserveForNextShell()) {
        // no next shell
        return false;
    }

    const auto depth = shells.size();
    const auto shell_start = neighbor_groups.size();

    for (auto& a : shells.back() | std::views::join) {
        a.expand(atoms);
        neighbor_groups.push_back(a.children);
    }
    shells.emplace_back(neighbor_groups.begin() + shell_start, neighbor_groups.end());
    return true;
}


/*
 * Sort substituents
 * State whether they are sorted
 * If sorted, is it "chiral" or "psuedochiral"
 * find unsorted subsets
 */

CMP sort_substituents(std::vector<CIPAtom>& substituents)
{
    std::sort(substituents.begin(), substituents.end(), cmp);
    return CMP::EQUAL;
}



// assume `Shells` understands how to build subsequent shell
// Currently templated on the function, could also template on an enum that
// controls what function is called
template <T cmp>
CMP rank_ligands(shells1, shells2)
{

    for (auto [shell1, shell2]: zip(shells1, shells2)) {
        if (shell1.empty() && !shell2.empty()) {
            return CMP::LESS;
        if (!shell1.empty() && shell2.empty()) {
            return CMP::GREATER;
        }
        // compare each set of substituents, starting with the neighbors
        // of the highest ranked neighbor in the previous group
        for (auto [substituents1, substituents2]: zip(shell1, shell2)) {
            if (substituents1.empty() && !substituents2.empty()) {
                return CMP::LESS;
            if (!substituents1.empty() && substituents2.empty()) {
                return CMP::GREATER;
            }
            // how many layers does this need to be?
            std::ranges::sort(substituents1, cmp);
            std::ranges::sort(substituents2, cmp);
            // I think it's the same ones as these
            res = cmp(substituents1, substituents2);
            if (res != CMP::EQUAL) {
                return res;
            }
        }
    }
}


CMP doCMP(auto this, auto that) {
    if (this < that) {
        return CMP::LESS;
    } else if (this == that) {
        return CMP::EQUAL;
    } else {
        return CMP::GREATER;
    }
}

// Rule 1a: Higher atomic number precedes lower.
CMP rule_1a(const CIPAtom& a1, const CIPAtom& a2) {
    return doCMP(a1.atomicNumber, a2.atomicNumber);
}

// Rule 1b (proposed): Lower root distance precedes higher root distance, where “root
// distance” is defined: (a) in the case of ring-closure duplicate nodes as the sphere of
// the duplicated atom; (b) in the case of multiple-bond duplicate nodes as the sphere of
// the atom to which the duplicate node is attached; and (c) in all other cases as the
// sphere of the atom itself.
bool rule_1b(const CIPAtom& a1, const CIPAtom& a2) {
    return doCMP(a1.depth, a2.depth);
}

bool rule_2(const CIPAtom& a1, const CIPAtom& a2) {
    return doCMP(a1.massNumber, a2.massNumber);
}

// Rule 3: When considering double bonds and planar tetraligand atoms, ‘seqcis’ = ‘Z’ precedes
// ‘seqtrans’ = ‘E’, and this precedes nonstereogenic double bonds.
bool rule_3(const CIPAtom& a1, const CIPAtom& a2) {
    // require that bond descriptors are set
    return doCMP(a1.bondDescriptor, a2.bondDescriptor);
}

// Rule 4a: Chiral stereogenic units precede pseudoasymmetric stereogenic units, and these precede
// nonstereogenic units.
bool rule_4a(const CIPAtom& a1, const CIPAtom& a2) {
    // require that atom descriptors are set
    return doCMP(a1.atomDescriptor / 2, a2.atomDescriptor / 2);
}


// Rule 4b: When two ligands have different descriptor pairs, then the one with the first chosen like
// descriptor pair has priority over the one with a corresponding unlike descriptor pair.
// this is the crazy one

// Rule 4c: ‘r’ precedes ‘s’ and ‘m’ precedes ‘p’.

// Rule 5: An atom or group with descriptor ‘R’, ‘M’, or ‘seqCis’ has priority over its enantiomorph ‘S’,
// ‘P’, or ‘seqTrans’.
//

// Rule 6 (proposed): An undifferentiated reference node has priority over any other
// undifferentiated node


// sort a range given a comparator. record whether it is
// fully sorted (or pseudo-sorted)
// simple bubble sort - we expect the ranges to be small (<4)
template <T cmp>
CMP shallow_cip_sort(auto& range) {
    using enum CMP;
    bool pseudo = false;
    bool unsortable = false;
    for (size_t i =0; i < range.size(); ++i) {
        for (size_t j= i + 1; j < range.size(); ++j) {
            auto res = cmp(range[i], range[j]);
            if (res > EQUAL) {
                std::swap(range[i], range[j]);
                // update child pointers here?
                for (auto& c: range[i].children) {
                    c.parent = &range[i];
                }
                for (auto& c: range[i].children) {
                    c.parent = &range[j];
                }
                updated = true;
            }
            if (res == EQUAL) {
                unsortable = true;
            }
            if (res == PSEUDO_LESS || res == PSEUDO_GREATER) {
                pseudo = true;
            }
        }
    }

    if (unsortable) {
        return EQUAL;
    } else if (pseudo) {
        return PSEUDO_LESS;
    } else {
        return LESS;
    }
}


template <auto... Rules, typename T>
CMP check_all_rules(const T& s1, const T& s2) {
    CMP res = CMP::EQUAL;

    // Fold expression: executes rank_ligands for each Rule in order
    // Short-circuits as soon as res != CMP::EQUAL
    ((res = rank_ligands<Rules>(s1, s2), res == CMP::EQUAL) && ...);

    return res;
}

CMP constitutional_rank_ligands(const Atom* src, const Atom* ligand1, const Atom* ligand2) {
    Shells shells1(src, ligand1);
    Shells shells2(src, ligand2);

    return check_all_rules<rule_1a, rule_1b, rule_2>(shells1, shells2);
}

// make this templated on the rule set
// take a max recursion depth parameter
Descriptor label(std::vector<std::pair<Atom*, Atom*> directors) {
    std::vector<Shells> ligands;
    for (auto [src, l]: directors) {
        ligands.emplace_back(src, l);
    }
    // set some flag that says a differentiator was pseudo
    std::ranges::sort(ligands, check_all_rules<rule_1a, rule_1b, rule_2>);
}
