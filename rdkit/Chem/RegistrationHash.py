"""
Generate a unique hash code for a molecule based on chemistry. If two
molecules are chemically "the same", they should have the same hash.

Used by Schrödinger's LiveDesign to determine if two molecules are
the same. LiveDesign makes changes to the molecule before molhash,
somewhat equivalent the steps available in
rdkit.Chem.MolStandardize.

Using molhash adds value beyond using SMILES because it:

* Ignores SMILES features that are not chemically meaningful
(e.g. atom map numbers)
* Canonicalizes enhanced stereochemistry groups. For example
`C[C@H](O)CC |&1:1|` and `C[C@@H](O)CC |&1:1|` have the same
molhash
* Canonicalizes S group data (for example, polymer data)

There are two hash schemes, the default, and one in which
tautomers are considered equivalent.


Copyright (C) 2022 Schrödinger, LLC

"""

import enum
import hashlib
import json
import logging
import re
from typing import Iterable
from typing import Optional

from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolHash

ATOM_PROP_MAP_NUMBER = 'molAtomMapNumber'

logger = logging.getLogger(__name__)

rdDepictor.SetPreferCoordGen(True)

EXTRA_ISOTOPE_REMOVAL_REGEX = re.compile(r'\[[1-3]0*([1-9]?[0-9]*[A-Z][a-z]?)@')
ENHANCED_STEREO_GROUP_REGEX = re.compile(r'((?:a|[&o]\d+):\d+(?:,\d+)*)')

ENHANCED_STEREO_GROUP_WEIGHTS = {
    Chem.StereoGroupType.STEREO_AND: 1000,
    Chem.StereoGroupType.STEREO_OR: 2000,
    Chem.StereoGroupType.STEREO_ABSOLUTE: 3000,
}

EMPTY_MOL_TAUTOMER_HASH = "_0_0"


class HashLayer(enum.Enum):
    """
    :cvar CANONICAL_SMILES: RDKit canonical SMILES (excluding enhanced stereo)
    :cvar ESCAPE: arbitrary other information to be incorporated
    :cvar FORMULA: a simple molecular formula for the molecule
    :cvar NO_STEREO_SMILES: RDKit canonical SMILES with all stereo removed
    :cvar NO_STEREO_TAUTOMER_HASH: the above tautomer hash lacking all stereo
    :cvar SGROUP_DATA: canonicalization of all SGroups data present
    :cvar TAUTOMER_HASH: SMILES-like representation for a generic tautomer form

    See SS-30145 for more documentation and example jupyter notebook
    """
    CANONICAL_SMILES = enum.auto()
    ESCAPE = enum.auto()
    FORMULA = enum.auto()
    NO_STEREO_SMILES = enum.auto()
    NO_STEREO_TAUTOMER_HASH = enum.auto()
    SGROUP_DATA = enum.auto()
    TAUTOMER_HASH = enum.auto()


@enum.unique
class HashScheme(enum.Enum):
    """
    Which hash layers to use to when deduplicating molecules

    Typically the "ALL_LAYERS" scheme is used, but some users may want
    the "TAUTOMER_INSENSITIVE_LAYERS" scheme.

    :cvar ALL_LAYERS: most strict hash scheme utilizing all layers
    :cvar STEREO_INSENSITIVE_LAYERS: excludes stereo sensitive layers
    :cvar TAUTOMER_INSENSITIVE_LAYERS: excludes tautomer sensitive layers
    """
    ALL_LAYERS = tuple(HashLayer)
    STEREO_INSENSITIVE_LAYERS = (
        HashLayer.ESCAPE,
        HashLayer.FORMULA,
        HashLayer.NO_STEREO_SMILES,
        HashLayer.NO_STEREO_TAUTOMER_HASH,
        HashLayer.SGROUP_DATA,
    )
    TAUTOMER_INSENSITIVE_LAYERS = (
        HashLayer.ESCAPE,
        HashLayer.FORMULA,
        HashLayer.NO_STEREO_TAUTOMER_HASH,
        HashLayer.SGROUP_DATA,
        HashLayer.TAUTOMER_HASH,
    )


def get_molhash(all_layers,
                hash_scheme: HashScheme = HashScheme.ALL_LAYERS) -> str:
    """
    Generate a molecular hash using a specified set of layers.

    :param mol: the molecule to generate the hash for
    :param hash_scheme: enum encoding information layers for the hash
    :return: hash for the given scheme constructed from the input layers
    """

    h = hashlib.sha1()
    for layer in hash_scheme.value:
        h.update(all_layers[layer].encode())
    return h.hexdigest()


def get_mol_layers(original_molecule: Chem.rdchem.Mol,
                   data_field_names: Optional[Iterable] = None,
                   escape: Optional[str] = None) -> set(HashLayer):
    """
    Generate layers of data about that could be used to identify a molecule

    :param original_molecule: molecule to obtain canonicalization layers from
    :param data_field_names: optional sequence of names of SGroup DAT fields which
       will be included in the hash.
    :param escape: optional field which can contain arbitrary information
    :return: dictionary of HashLayer enum to calculated hash
    """
    # Work on a copy with all non-stereogenic Hydrogens removed
    mol = _remove_unnecessary_hs(original_molecule,
                                 preserve_stereogenic_hs=True)
    strip_atom_map_labels(mol)

    formula = rdMolHash.MolHash(mol, rdMolHash.HashFunction.MolFormula)
    cxsmiles, canonical_mol = canonicalize_stereo_groups(mol)
    tautomer_hash = get_stereo_tautomer_hash(canonical_mol)

    canonical_smiles = get_canonical_smiles(cxsmiles)
    sgroup_data = canonicalize_sgroups(canonical_mol,
                                       dataFieldNames=data_field_names)

    no_stereo_tautomer_hash, no_stereo_smiles = get_no_stereo_layers(
        canonical_mol)

    return {
        HashLayer.CANONICAL_SMILES: canonical_smiles,
        HashLayer.ESCAPE: escape or "",
        HashLayer.FORMULA: formula,
        HashLayer.NO_STEREO_SMILES: no_stereo_smiles,
        HashLayer.NO_STEREO_TAUTOMER_HASH: no_stereo_tautomer_hash,
        HashLayer.SGROUP_DATA: sgroup_data,
        HashLayer.TAUTOMER_HASH: tautomer_hash,
    }


def strip_atom_map_labels(mol):
    for at in mol.GetAtoms():
        if at.HasProp(ATOM_PROP_MAP_NUMBER):
            at.ClearProp(ATOM_PROP_MAP_NUMBER)


def _remove_unnecessary_hs(rdk_mol, preserve_stereogenic_hs=False):
    """
    removes hydrogens that are not necessary for the registration hash, and
    preserves hydrogen isotopes
    """
    remove_hs_params = Chem.RemoveHsParameters()
    remove_hs_params.updateExplicitCount = True
    remove_hs_params.removeDefiningBondStereo = not preserve_stereogenic_hs
    edited_mol = Chem.rdmolops.RemoveHs(rdk_mol,
                                        remove_hs_params,
                                        sanitize=False)
    edited_mol.UpdatePropertyCache(False)
    return edited_mol


def get_stereo_tautomer_hash(molecule):
    if molecule.GetNumAtoms() == 0:
        return EMPTY_MOL_TAUTOMER_HASH

    # SHARED-7909: workaround for https://github.com/rdkit/rdkit/issues/4234
    # This can be removed once we update to an RDKit version without this bug.
    no_h_mol = _remove_unnecessary_hs(molecule)
    no_h_mol.UpdatePropertyCache(False)

    # setting useCxSmiles param value to always include enhanced stereo info
    useCxSmiles = True
    hash_with_cxExtensions = rdMolHash.MolHash(
        no_h_mol, rdMolHash.HashFunction.HetAtomTautomer, useCxSmiles)

    # since the cxSmiles can include anything, we want to only include
    # enhanced stereo information
    canonical_smiles = get_canonical_smiles(hash_with_cxExtensions)
    return canonical_smiles


def get_canonical_smiles(cxsmiles):
    smiles_parts = (p.strip() for p in cxsmiles.split('|'))
    smiles_parts = [p for p in smiles_parts if p]
    if not smiles_parts:
        return '', ''
    elif len(smiles_parts) > 2:
        raise ValueError('Unexpected number of fragments in canonical CXSMILES')

    canonical_smiles = smiles_parts[0]
    stereo_groups = ''
    if len(smiles_parts) == 2:
        # note: as with many regex-based things, this is fragile
        groups = ENHANCED_STEREO_GROUP_REGEX.findall(smiles_parts[1])
        if groups:
            # We might have other CXSMILES extensions that aren't stereo groups
            stereo_groups = f"|{','.join(sorted(groups))}|"

    if stereo_groups:
        return canonical_smiles + " " + stereo_groups

    return canonical_smiles


def get_no_stereo_layers(mol):
    # SHARED-7909: Strip all Hs, including stereogenic ones.
    no_stereo_mol = _remove_unnecessary_hs(mol)
    no_stereo_mol.UpdatePropertyCache(False)
    Chem.rdmolops.RemoveStereochemistry(no_stereo_mol)
    no_stereo_tautomer_hash = rdMolHash.MolHash(
        no_stereo_mol, rdMolHash.HashFunction.HetAtomTautomer)
    no_stereo_smiles = rdMolHash.MolHash(no_stereo_mol,
                                         rdMolHash.HashFunction.CanonicalSmiles)
    return (no_stereo_tautomer_hash, no_stereo_smiles)


def get_canonical_atom_ranks_and_bonds(mol, useSmilesOrdering=True):
    """
    returns a 2-tuple with:

    1. the canonical ranks of a molecule's atoms
    2. the bonds expressed as (canonical_atom_rank_1,canonical_atom_rank_2) where
       canonical_atom_rank_1 < canonical_atom_rank_2

    If useSmilesOrdering is True then the atom indices here correspond to the order of
    the atoms in the canonical SMILES, otherwise just the canonical atom order is used.
    useSmilesOrdering=True is a bit slower, but it allows the output to be linked to the
    canonical SMILES, which can be useful.

    """
    if mol.GetNumAtoms() == 0:
        return [], []
    if not useSmilesOrdering:
        atRanks = list(Chem.CanonicalRankAtoms(mol))
    else:
        smi = Chem.MolToSmiles(mol)
        ordertxt = mol.GetProp("_smilesAtomOutputOrder")
        smiOrder = [int(x) for x in ordertxt[1:-1].split(",") if x]
        atRanks = [0] * len(smiOrder)
        for i, idx in enumerate(smiOrder):
            atRanks[idx] = i
    # get the bonds in canonical form
    bndOrder = []
    for bnd in mol.GetBonds():
        bo = atRanks[bnd.GetBeginAtomIdx()]
        eo = atRanks[bnd.GetEndAtomIdx()]
        if bo > eo:
            bo, eo = eo, bo
        bndOrder.append((bo, eo))
    return atRanks, bndOrder


def canonicalize_data_sgroup(sg,
                             atRanks,
                             bndOrder,
                             fieldNames=None,
                             sortAtomOrder=True):
    """
    NOTES: if sortAtomOrder is true then the atom list will be sorted.
    This assumes that the order of the atoms in that list is not important

    """
    if sg.GetProp("TYPE") != "DAT" or not sg.HasProp("FIELDNAME"):
        return None

    canonicalizableFieldNames = fieldNames or ["Atrop"]
    fieldName = sg.GetProp("FIELDNAME")
    if fieldName not in canonicalizableFieldNames:
        return None

    data = sg.GetStringVectProp("DATAFIELDS")
    if len(data) > 1:
        raise ValueError(
            "cannot canonicalize data groups with multiple data fields")
    data = data[0]
    ats = tuple(atRanks[x] for x in sg.GetAtoms())
    if sortAtomOrder:
        ats = tuple(sorted(ats))
    bnds = tuple(bndOrder[x] for x in sg.GetBonds())

    res = dict(fieldName=fieldName, atom=ats, bonds=bnds, value=data)
    return res


def getCanonicalBondRep(bond, atomRanks):
    aid1 = bond.GetBeginAtomIdx()
    aid2 = bond.GetEndAtomIdx()
    if atomRanks[aid1] > atomRanks[aid2] or (atomRanks[aid1] == atomRanks[aid2]
                                             and aid1 > aid2):
        aid1, aid2 = aid2, aid1
    return (aid1, aid2)


def canonicalize_sru_sgroup(mol, sg, atRanks, bndOrder, sortAtomAndBondOrder):
    """
    NOTES: if sortAtomAndBondOrder is true then the atom and bond lists will be sorted.
    This assumes that the ordering of those lists is not important

    """
    if sg.GetProp("TYPE") != "SRU":
        return None
    ats = tuple(atRanks[x] for x in sg.GetAtoms())
    bnds = tuple(bndOrder[x] for x in sg.GetBonds())
    if sortAtomAndBondOrder:
        ats = tuple(sorted(ats))
        bnds = tuple(sorted(bnds))
    props = sg.GetPropsAsDict()
    res = dict(
        type="SRU",
        atoms=ats,
        bonds=bnds,
        index=props.get("index", 0),
        connect=props.get("CONNECT", ""),
        label=props.get("LABEL", ""),
    )
    if "PARENT" in props:
        res["PARENT"] = props["PARENT"]
    if "XBHEAD" in props:
        xbhbonds = tuple(
            getCanonicalBondRep(mol.GetBondWithIdx(x), atRanks)
            for x in props["XBHEAD"])
        if sortAtomAndBondOrder:
            xbhbonds = tuple(sorted(xbhbonds))
        res["XBHEAD"] = xbhbonds
    if "XBCORR" in props:
        xbcorrbonds = tuple(
            getCanonicalBondRep(mol.GetBondWithIdx(x), atRanks)
            for x in props["XBCORR"])
        if len(xbcorrbonds) % 2:
            raise ValueError("XBCORR should have 2N bonds")
        if sortAtomAndBondOrder:
            # these are pairs, so we need to sort them as such
            tmp = []
            for i in range(0, len(xbcorrbonds), 2):
                b1, b2 = xbcorrbonds[i], xbcorrbonds[i + 1]
                if b1 > b2:
                    b1, b2 = b2, b1
                tmp.append((b1, b2))
            xbcorrbonds = tuple(sorted(tmp))
        res["XBCORR"] = xbcorrbonds
    return res


def canonicalize_cop_sgroup(sg, atRanks, sortAtomAndBondOrder):
    """
    NOTES: if sortAtomAndBondOrder is true then the atom and bond lists will be sorted.
    This assumes that the ordering of those lists is not important
    """
    if sg.GetProp("TYPE") != "COP":
        return None
    ats = tuple(atRanks[x] for x in sg.GetAtoms())
    if sortAtomAndBondOrder:
        ats = tuple(sorted(ats))
    props = sg.GetPropsAsDict()
    # we are intentionally not doing the label here. Marvin adds one, biovia draw does not
    res = dict(type="COP", atoms=ats, index=props.get("index", 0))
    return res


def canonicalize_sgroups(mol, dataFieldNames=None, sortAtomAndBondOrder=True):
    """
    NOTES: if sortAtomAndBondOrder is true then the atom and bond lists will be sorted.
    This assumes that the ordering of those lists is not important
    """
    dataFieldNames = dataFieldNames or ["Atrop"]
    atRanks, bndOrder = get_canonical_atom_ranks_and_bonds(mol)
    res = []
    for sg in Chem.GetMolSubstanceGroups(mol):
        lres = None
        if sg.GetProp("TYPE") == "DAT":
            lres = canonicalize_data_sgroup(sg, atRanks, bndOrder,
                                            dataFieldNames,
                                            sortAtomAndBondOrder)
        elif sg.GetProp("TYPE") == "SRU":
            lres = canonicalize_sru_sgroup(mol, sg, atRanks, bndOrder,
                                           sortAtomAndBondOrder)
        elif sg.GetProp("TYPE") == "COP":
            lres = canonicalize_cop_sgroup(sg, atRanks, sortAtomAndBondOrder)

        if lres is not None:
            res.append(lres)
    if len(res) > 1:
        # we need to sort them, but dicts don't let you do that directly:
        tres = sorted(tuple(x.items()) for x in res)
        res = tuple(dict(x) for x in tres)
        # update "index" and "PARENT" props
        idxmap = {}
        for i, itm in enumerate(res):
            if "index" in itm:
                idxmap[itm["index"]] = i + 1
                itm["index"] = i + 1
        for i, itm in enumerate(res):
            if "PARENT" in itm:
                itm["PARENT"] = idxmap[itm["PARENT"]]

    return json.dumps(res)


class EnhancedStereoUpdateMode(enum.Enum):
    ADD_WEIGHTS = enum.auto()
    REMOVE_WEIGHTS = enum.auto()


def update_enhanced_stereo_group_weights(mol, mode):
    if mode == EnhancedStereoUpdateMode.ADD_WEIGHTS:
        factor = 1
    elif mode == EnhancedStereoUpdateMode.REMOVE_WEIGHTS:
        factor = -1
    else:
        raise ValueError('Invalid Enhanced Stereo weight update mode')

    isotopesModified = False
    stgs = mol.GetStereoGroups()
    for stg in stgs:
        stgt = stg.GetGroupType()
        weight = factor * ENHANCED_STEREO_GROUP_WEIGHTS[stgt]
        for at in stg.GetAtoms():

            # Make sure the isotope is reasonable and is not present in more than
            # one stereo group, and we are not messing it up
            isotope = at.GetIsotope()
            if mode == EnhancedStereoUpdateMode.ADD_WEIGHTS and isotope >= 100:
                raise ValueError(
                    f'Atom {at.GetIdx()} is an unusual isotope ({isotope})')

            at.SetIsotope(isotope + weight)
            isotopesModified = True

    return mol, isotopesModified


def canonicalize_stereo_groups(mol):
    """
    Returns canonical CXSmiles and the corresponding molecule with the
    stereo groups canonicalized.

    The RDKit canonicalization code does not currently take stereo groups into
    account. We work around that by using EnumerateStereoisomers() to generate
    all possible instances of the molecule's stereogroups and then lexically
    compare the CXSMILES of those.
    """

    if not len(mol.GetStereoGroups()):
        return Chem.MolToCXSmiles(mol), mol

    mol = Chem.Mol(mol)

    # add ring info if not initialized yet
    Chem.FastFindRings(mol)

    # To solve the problem we add isotope tags to the atoms involved in stereo
    # groups. These allow the canonicalization to tell the difference between
    # AND and OR (and have the happy side effect of allowing the
    # presence/absence of a stereo group to affect the canonicalization)
    mol, isotopesModified = update_enhanced_stereo_group_weights(
        mol, EnhancedStereoUpdateMode.ADD_WEIGHTS)

    # We're going to be generating canonical SMILES here anyway, so skip the
    # "unique" option for the sake of efficiency
    opts = EnumerateStereoisomers.StereoEnumerationOptions()
    opts.onlyStereoGroups = True
    opts.unique = False

    res = None
    res_csmi = ''
    for isomer in EnumerateStereoisomers.EnumerateStereoisomers(mol, opts):
        # We need to generate the canonical CXSMILES for the molecule with
        # the isotope tags.
        csmi = Chem.MolToCXSmiles(isomer)
        if res is None or csmi < res_csmi:
            res = isomer
            res_csmi = csmi

    res_csmi = EXTRA_ISOTOPE_REMOVAL_REGEX.sub(r'[\1@', res_csmi)

    if isotopesModified:
        res, _ = update_enhanced_stereo_group_weights(
            res, EnhancedStereoUpdateMode.REMOVE_WEIGHTS)

    return res_csmi, res