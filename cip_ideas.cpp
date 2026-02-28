
#include <QtCore/qnamespace.h>
#include <algorithm>
#include <ranges>
#include <iterator>

struct Mol {

vector<Atom> atoms;

vector<int> original_indices;

CIPMol(ROMol& m) {
    // populate atoms
    // reorder atoms by connectivity
}

};



void label(ROMol mol) {
    CIPMol cip_mol(mol);

}

struct Atom {
int atomic_number;
int isotope;
};

auto span(v, indices) {
    return std::span(v.begin() + indices[0], v.begin() + indices[1]);
}

void sort(v, indices, cmp)
{
    std::sort(v.begin() + indices[0], v.begin() + indices[1], cmp);
}


template <T cmp>
rank_ligands(mol, src, ligand1, ligand2)
{
    // These should be bonds because we only go one way
    vector frontier1 = {{src, ligand1}};
    vector frontier_groups1 = {{0, 1}};
    vector frontier2 = {{src, ligand2}};
    vector frontier_groups2 = {{0, 1}};

    // two sets of frontiers because we always have "current"
    // and "next". Swap 'em at each iteration
    decltype(frontier1) next_frontier1;
    decltype(frontier_groups1) next_frontier_groups1;
    decltype(frontier2) next_frontier2;
    decltype(frontier_groups2) next_frontier_groups2;

    constexpr MAXDEPTH = 12;

    auto res = SAME;
    auto depth = 0;
    while (res == SAME && depth < MAXDEPTH) {
        for (auto s1, s2: zip(frontier1, frontier2)) {
            auto a1 = span(this, s1);
            auto a2 = span(this, s2);
            std::ranges::sort(a1, cmp);
            std::ranges::sort(a2, cmp);
            res = cmp(a1, a2);
            if (res != 0) {
                return res;
            }
            // expand the shells
            for (auto a: a1) {
                auto start = next_frontier1.size();
                for (auto n: a.neighbors) {
                    next_frontier1.push_back(n);
                    // Add a depth feature
                    // check whether the node is a "duplicate" or "real" node
                    // check here about "virtualization"
                }
                frontier_groups1.push_back({start, next_frontier1.size()});
            }
        }

        // swap our working buffer for the next shell, clear space for the subsequent
        // shell
        std::swap(frontier1, next_frontier1);
        std::swap(frontier2, next_frontier2);
        std::swap(frontier_groups1, next_frontier_groups1);
        std::swap(frontier_groups2, next_frontier_groups2);
        next_frontier1.resize(0);
        next_frontier2.resize(0);
        next_frontier_groups1.resize(0);
        next_frontier_groups2.resize(0);
    }


}


enum class CMP {
LESS, PSEUDO_LESS, EQUAL, PSEUDO_GREATER, GREATER
}
static constexpr auto EMPTY

struct CIPAtom {
uchar atomic_number; // this may be fractional due to averaging
double mass_number; // this may be fractional due to averaging
size_t depth;
Descriptor descriptor;
char flags; // ring duplicate, bond duplicate, implicit hydrogen

size_t parent = EMPTY;
size_t idx = EMPTY;
std::span<CIPAtom> children;
Shells& mol;

// link to the real atom

bool isVisited(size_t i) const {
    auto a = this;
    while (a.idx != i && a.parent != EMPTY) {
        a = &(a.mol.atoms[a.parent]);
    }
    return a.idx == i;
}
};

// actual atoms to iterate over.
// maybe should be a span of Atom?
stuct NeighborList
{

};

stuct NextShell
{
bool empty() const;

// some range adapter work here
NeighborList* begin();
NeighborList* end();
};

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
    Shells(mol, src, ligand) {

    }

    // iterate over the shells out from something
    ShellIterator* begin() {return ShellIterator(this, false);}
    ShellIterator* end() {return ShellIterator(this, true);}

    std::span<std::span<CIPAtom>>* getShell(size_t idx);
    // allow re-use of the memory
    // assign
    // clear

private:
    // populate the next shell and fix spans
    bool reserveForNextShell();
    bool makeNextShell();

    std::vector<char> ever_visited;

    std::vector<CIPAtom> atoms; // all atoms
    std::vector<std::span<CIPAtom>> neighbor_groups; // groups of neighbors of some atom
    std::vector<std::array<size_t, 2> neighbor_groups_idx; // index-based for growing
    std::vector<std::span<std::span<CIPAtom>>> shells; // iterable of neighbor groups
    std::vector<std::array<size_t, 2>> shells_idx; // index-based for growing

};

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
    size_t new_neighbor_group_count = 0;
    for (auto& a : shells.back() | std::views::join) {
        neighbors = a.neighbors.size() - 1; // ignore parent
        new_neighbor_group_count += 1 ? neighbors != 0 : 1;
        new_atom_count += neighbors;
    }
    if (new_atom_count == 0) {
        return false;
    }
    const auto atoms_ptr = std::static_cast<void*>(atoms.data());
    atoms.reserve(atoms.size() + new_atom_count);

    const auto neighbor_groups_ptr = std::static_cast<void*>(neighbor_groups.data());
    neighbor_groups.reserve(neighbor_groups.size() + new_neighbor_group_count);
    neighbor_groups_idx.reserve(neighbor_groups_idx.size() + new_neighbor_group_count);
    // reallocation, need to update span pointers
    if (atoms_ptr != std::static_cast<void*>(atoms.data())) {
        for (size_t i =0; i < neighbor_groups.size(); ++i) {
            auto [b, e] = neighbor_groups_idx[i];
            neighbor_groups[i] = std::span(atoms.begin() +b, atoms.begin() +e);
            // reset children of the parent
            if (!neighbor_groups.empty()) {
                auto& parent = atoms[neighbor_groups[i].front().parent];
                parent.children = neighbor_groups[i];
            }
        }
    }

    shells.reserve(shells.size() + new_neighbor_group_count);
    shells_idx.reserve(shells_idx.size() + new_neighbor_group_count);
    // reallocation, need to update span pointers
    if (neighbor_groups_ptr != std::static_cast<void*>(neighbor_groups.data())) {
        for (size_t i =0; i < shells.size(); ++i) {
            auto [b, e] = shells_idx[i];
            shells[i] = std::span(neighbor_groups.begin() +b, neighbor_groups.begin() +e);
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
        const auto next_group_start = atoms.size();

        for (auto n: a.neighbors) {
            if (n == a.parent) {continue;}
            // check for loop closure and other fake nodes
            if (ever_visited[n.atom->getIdx()]) {
                // check the path back
            }
            //
            atoms.emplace_back {
                n
                .depth=depth
            };
        }
        // add implicit H nodes
        neighbor_groups_idx.emplace_back(next_group_start, atoms.size());
        neighbor_groups.emplace_back(atoms.begin() + next_group_start, atoms.end());
    }
    shells_idx.emplace_back(shell_start, neighbor_groups.size() - shell_start);
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
//
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
            std::ranges::sort(substituents1, cmp);
            std::ranges::sort(substituents2, cmp);
            res = cmp(substituents1, substituents2);
            if (res != CMP::EQUAL) {
                return res;
            }
        }
    }
}


bool rule_1a(Atom a1, Atom a2) {
    return a1.atomic_number < a2.atomic_number;
}
bool rule_1b(Atom a1, Atom a2) {
    if (a1.atomic_number < a2.atomic_number) {
        return true;
    } else if (a1.atomic_number > a2.atomic_number) {
        return false;
    }
    return a1.mass_number < a2.mass_number;
}


template <auto... Rules, typename T>
CMP check_all_rules(const T& s1, const T& s2) {
    CMP res = CMP::EQUAL;

    // Fold expression: executes rank_ligands for each Rule in order
    // Short-circuits as soon as res != CMP::EQUAL
    ((res = rank_ligands<Rules>(s1, s2), res == CMP::EQUAL) && ...);

    return res;
}

CMP rank_ligands(mol, src, ligand1, ligand2) {
    Shells shells1(mol, src, ligand1);
    Shells shells2(mol, src, ligand2);

    return check_all_rules<rule_1a, rule_1b, rule_2>(shells1, shells2);
}
