/*
 * vector<CIPAtom> atoms
 * children should be a span
 * comparing across the groups of children in out-most depth
 *  decides order.
 * only compare ties from earlier generations
 * need to be able to track parity stacks
 * tie resolution re-sorts children
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
    // do the constitutional rankings
    // do the remaining stereo bonds
    // - add auxiliary descriptors, out to in
    // - label atoms
    // do the remaining stereo atoms
    // - add auxiliary descriptors, out to in
    // - label atoms
    //

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

bool is_pseudo(ATOMDESCRIPTOR descriptor) {
    using enum ATOMDESCRIPTOR;
    return descriptor == r || descriptor == s;
}
// bool is_pseudo(BONDDESCRIPTOR descriptor) {
//     using enum BONDDESCRIPTOR;
//     return descriptor == m || descriptor == p;
// }
bool is_real(ATOMDESCRIPTOR descriptor) {
    using enum ATOMDESCRIPTOR;
    return descriptor == R || descriptor == S;
}
bool is_real(BONDDESCRIPTOR descriptor) {
    using enum BONDDESCRIPTOR;
    return descriptor == Z || descriptor == E;
}


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
                // add bond stereo info
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
    // Make a CIPAtom reference to another CIPAtom in an existing, fully expanded graph.
    // This happens the second (and Nth) expansion.
    //
    // Probably want a limit here so that it doesn't recurse indefinitely.
}
CIPAtom::CIPAtom(CIPAtom* parent, Atom* reference, bool bondDuplicate):
    depth(parent->depth + 1), realAtom(reference), parent(parent), bondDuplicate(bondDuplicate)
{
    // Make a CIPAtom based on a regular RDKit atom (or implicit H).
    // This happens the first epansion of at shell tree.
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
    std::vector<std::array<CIPAtom*, 2>> ties;
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

    for (auto& [l, r]: ties) {
        l += new_atoms_offset;
        r += new_atoms_offset;
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
CMP rule_1b(const CIPAtom& a1, const CIPAtom& a2) {
    return doCMP(a1.depth, a2.depth);
}

CMP rule_2(const CIPAtom& a1, const CIPAtom& a2) {
    return doCMP(a1.massNumber, a2.massNumber);
}

// Rule 3: When considering double bonds and planar tetraligand atoms, ‘seqcis’ = ‘Z’ precedes
// ‘seqtrans’ = ‘E’, and this precedes nonstereogenic double bonds.
CMP rule_3(const CIPAtom& a1, const CIPAtom& a2) {
    // require that bond descriptors are set
    return doCMP(a1.bondDescriptor, a2.bondDescriptor);
}

// Rule 4a: Chiral stereogenic units precede pseudoasymmetric stereogenic units, and these precede
// nonstereogenic units.
CMP rule_4a(const CIPAtom& a1, const CIPAtom& a2) {
    // require that atom descriptors are set
    // should this also do bond descriptors? I think so
    return doCMP(a1.atomDescriptor / 2, a2.atomDescriptor / 2);
}


// Rule 4b: When two ligands have different descriptor pairs, then the one with the first chosen like
// descriptor pair has priority over the one with a corresponding unlike descriptor pair.
// this is the crazy one

// Rule 4c: ‘r’ precedes ‘s’ and ‘m’ precedes ‘p’.
CMP rule_4c(const CIPAtom& a1, const CIPAtom& a2) {
    if (is_pseudo(a1.atomDescriptor) && is_pseudo(a2.atomDescriptor)) {
        return doCMP(a1.atomDescriptor, a2.atomDescriptor);
    }
    // if (is_pseudo(a1.bondDescriptor) && is_pseudo(a2.bondDescriptor)) {
    //     return doCMP(a1.bondDescriptor, a2.bondDescriptor);
    // }
    return CMP::EQUAL;
}

// Rule 5: An atom or group with descriptor ‘R’, ‘M’, or ‘seqCis’ has priority over its enantiomorph ‘S’,
// ‘P’, or ‘seqTrans’.
CMP rule_5(const CIPAtom& a1, const CIPAtom& a2) {
    if (is_real(a1.atomDescriptor) && is_real(a2.atomDescriptor)) {
        return doCMP(a1.atomDescriptor, a2.atomDescriptor);
    }
    if (is_real(a1.bondDescriptor) && is_real(a2.bondDescriptor)) {
        return doCMP(a1.bondDescriptor, a2.bondDescriptor);
    }
    return CMP::EQUAL;
}

// Rule 6 (proposed): An undifferentiated reference node has priority over any other
// undifferentiated node
CMP rule_6(const CIPAtom& a1, const CIPAtom& a2) {
    // is this, like, has a copy somewhere?
    return CMP::EQUAL;
}


// sort a range given a comparator. record whether it is
// fully sorted (or pseudo-sorted)
// simple bubble sort - we expect the ranges to be small (<4)
template <T cmp>
CMP shallow_cip_sort(auto& range, auto& ties) {
    using enum CMP;
    bool pseudo = false;
    bool unsortable = false;
    vector<array<int, 2> tt;
    for (size_t i =0; i < range.size(); ++i) {
        for (size_t j= i + 1; j < range.size(); ++j) {
            auto res = cmp(range[i], range[j]);
            if (res > EQUAL) {
                swap_em(range[i], range[j]);
            }
            if (res == EQUAL) {
                tt.emplace_back(range[i].idx, range[j].idx);
                unsortable = true;
            }
            if (res == PSEUDO_LESS || res == PSEUDO_GREATER) {
                pseudo = true;
            }
        }
    }

    if (unsortable) {
        for (auto [l, r]: tt) {
            auto a = std::ranges::find(range, l);
            auto b = std::ranges::find(range, r);
            ties.emplace_back(a, b);
        }
        return EQUAL;
    } else if (pseudo) {
        return PSEUDO_LESS;
    } else {
        return LESS;
    }
}

auto children(CIPAtom* source)
{
    std::vector q{source->children};
    while (!q.empty()) {
        auto r = q.pop_front();
        co_yield r;
        for (const auto& a: r) {
            q.push_back(r.children);
        }
    }
}

auto child_groups_at_depth(CIPAtom* source, size_t depth)
{
    std::vector q{source};
    while (!q.empty()) {
        auto r = q.pop_front();
        if (!r->expanded) {
            continue;
        }
        if (r && r->depth + 1 == depth) {
            co_yield r->children;
        }
        else {
            for (const auto& c: r->children) {
                q.push_back(&c);
            }
        }
    }
}

auto zip_longest(auto& a, auto& b) {
    const auto short_l = std::min(a.size(), b.size());
    for (size_t i=0; i < short_l; ++i) {
        co_yield {&a[i], &b[i]};
    }
    if (a.size() > b.size()) {
        for (auto i=short_l, i < a.size(); ++i) {
            co_yield {&a[i], nullptr};
        }
    } else if (a.size() < b.size()) {
        for (auto i=short_l, i < b.size(); ++i) {
            co_yield {nullptr, &b[i]};
        }
    }
}

void swap_em(CIPAtom& a, CIPAtom& b) {
    std::swap(a, b);
    // update child pointers
    for (auto& c: a.children) {
        c.parent = &a;
    }
    for (auto& c: b.children) {
        c.parent = &b;
    }
}


// resolve ties within a group
// changes order within that group
template <T cmp>
CMP resolve_tie(std::array<CIPAtom*, 2>& tie, size_t depth)
{
    for (const auto& [group1, group2] : zip_longest(child_groups_at_depth(tie[0], depth), child_groups_at_depth(tie[1], depth))) {
        if (!group1 && !group2) {
            continue;
        }
        if ((!group1 ||
            group1.empty()) && (group2 && !group2->empty())) {
            return CMP::LESS;
        if (!group1.empty() && group2.empty()) {
            return CMP::GREATER;
        }
        auto res = cmp(a1, a2);
        if (res != CMP::EQUAL) {
            if (res < CMP::EQUAL) {
                swap_em(*tie[0], *tie[0]);
            }
            return res;
        }
    }
    return CMP::EQUAL;
}

template <T cmp>
CMP deep_cip_sort(auto& shells) {
    // try to sort
    // if a sort fails, expand to the next shell
    // once the next shell is sorted, break ties up the ladder?

    for (auto& shell: shells) {
        if (shell1.empty() && !shell2.empty()) {
            return CMP::LESS;
        if (!shell1.empty() && shell2.empty()) {
            return CMP::GREATER;
        }
        for (auto& group: shell) {
            // make the ties here
            shallow_cip_sort<cmp>(group, shells.ties);
        }
        CMP resolution;
        for (auto& tie: std::ranges::reverse(shells.ties)) {
            resolution = resolve_tie(tie);
            if (resolution != CMP::EQUAL) {
                shells.ties.remove(tie);
            }
        }
        if (resolution != CMP::EQUAL) {
            return resolution;
        }
    }
    return CMP::EQUAL;
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
