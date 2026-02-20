


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

{
    
    
    auto vector<atom> this;
    auto vector<std::span<atom>> this;
    auto vector<atom> that;
    auto vector<std::span<atom>> that;
    for (auto s0, s1: zip(this, that)) {
        auto a0 = span(this, s0);
        auto a1 = span(this, s1);
        std::ranges::sort(a0, cmp);
        std::ranges::sort(a1, cmp);
        res = cmp(a0, a1);
        if (res != 0) {
            return res;
        }
        // expand the shells
    }
    
}
