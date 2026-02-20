


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
