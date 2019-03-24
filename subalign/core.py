#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import chimera
import cStringIO as StringIO
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform
try:
    from rdkit import Chem
    from rdkit.Chem import MolFromPDBBlock, FastFindRings, SanitizeMol, Mol, RWMol, Atom, Conformer
    from rdkit.Chem.rdMolAlign import GetAlignmentTransform, GetO3A, AlignMol
    from rdkit.Chem.AllChem import MMFFGetMoleculeProperties, GetBestRMS
except ImportError:
    raise chimera.UserError("RDKit must be installed to use this extension.")

import numpy as np

# The RDKit_BondType dictionary is defined to convert float bond orders into RDKit bond types
RDKit_BondType = {
    1.0 : Chem.BondType.SINGLE,
    2.0 : Chem.BondType.DOUBLE,
    3.0 : Chem.BondType.TRIPLE,
    4.0 : Chem.BondType.QUADRUPLE,
    5.0 : Chem.BondType.QUINTUPLE,
    6.0 : Chem.BondType.HEXTUPLE,
    1.5 : Chem.BondType.ONEANDAHALF,
    2.5 : Chem.BondType.TWOANDAHALF,
    3.5 : Chem.BondType.THREEANDAHALF,
    4.5 : Chem.BondType.FOURANDAHALF,
    5.5 : Chem.BondType.FIVEANDAHALF
}

def _chimera_to_rdkit(molecule, sanitize=True):
    io = StringIO.StringIO()
    xform = molecule.openState.xform
    atoms = [a for a in molecule.atoms if not a.element.isMetal and a.element.number > 1]
    if not atoms:
        raise chimera.UserError("Molecule does not contain meaningful atoms!")
    molecule_copy = molecule_from_atoms(molecule, atoms)
    chimera.pdbWrite([molecule_copy], xform, io)
    io.seek(0)
    rdkit_mol = MolFromPDBBlock(io.getvalue(), False, sanitize)
    io.close()
    molecule_copy.destroy()
    return rdkit_mol


#Function to fix Nitro groups from chimera molecules without explicit formal charges
def _fix_nitro(molecule):
    """
    molecule : rdkit molecule
    """
    for atom in molecule.GetAtoms():
        valence = 0.0
        for bond in atom.GetBonds():
            valence += bond.GetBondTypeAsDouble()

        if atom.GetAtomicNum() == 7:
            if atom.GetFormalCharge():
                continue
            aromHolder = atom.GetIsAromatic()
            atom.SetIsAromatic(0)
            if valence == 4:
                for neighbor in atom.GetNeighbors():
                    if (neighbor.GetAtomicNum() == 8) and (molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE):
                        atom.SetFormalCharge(1)
                        neighbor.SetFormalCharge(-1)
                        break
                atom.SetIsAromatic(aromHolder)


def _transform_molecule(molecule, xform):
    new_xform = molecule.openState.xform
    new_xform.premultiply(xform)
    molecule.openState.xform = new_xform


def untransformed_rmsd(reference, probe, sanitize=True, uniquify=False,
                       reflect=True, method='sub', **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    if method in ('sub', 'best'):
        matches = rdk_reference.GetSubstructMatches(rdk_probe, uniquify=uniquify)
    elif method == 'o3a':
        rdk_reference2 = _chimera_to_rdkit(reference, sanitize=sanitize)
        rdk_probe2 = _chimera_to_rdkit(probe, sanitize=sanitize)
        FastFindRings(rdk_reference)
        FastFindRings(rdk_probe)
        reference_params = MMFFGetMoleculeProperties(rdk_reference)
        probe_params = MMFFGetMoleculeProperties(rdk_probe)
        o3a = GetO3A(rdk_probe2, rdk_reference2, probe_params, reference_params,
                     maxIters=0, reflect=reflect)
        matches = o3a.Matches()
    else:
        raise chimera.UserError('`method` must be sub, best or o3a')
    if not matches:
        raise chimera.UserError('Could not find any matches.')
    maps = [list(enumerate(match)) for match in matches]
    best_rmsd = 1000.
    for atom_map in maps:
        rmsd = AlignMol(rdk_probe, rdk_reference, -1, -1, atomMap=atom_map,
                       maxIters=0, reflect=reflect)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
    return best_rmsd


def align(reference, probe, transform=True, sanitize=True, maxIters=50, reflect=False, **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    rmsd, xform = GetAlignmentTransform(rdk_probe, rdk_reference,
                                        maxIters=maxIters, reflect=reflect)
    if transform:
        _transform_molecule(probe, chimera_xform(xform[:3]))
    return rmsd


def align_o3a(reference, probe, transform=True, sanitize=True,maxIters=50, reflect=False, **kwargs):

    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    FastFindRings(rdk_reference)
    FastFindRings(rdk_probe)
    reference_params = MMFFGetMoleculeProperties(rdk_reference)
    probe_params = MMFFGetMoleculeProperties(rdk_probe)
    o3a = GetO3A(rdk_probe, rdk_reference, prbPyMMFFMolProperties=probe_params,
                 refPyMMFFMolProperties=reference_params, maxIters=maxIters,
                 reflect=reflect)
    rmsd, xform = o3a.Trans()
    if transform:
        _transform_molecule(probe, chimera_xform(xform[:3]))
    return rmsd


def align_best(reference, probe, transform=True, sanitize=True, ignore_warnings=False,
               maxIters=50, reflect=False, **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    matches = rdk_reference.GetSubstructMatches(rdk_probe, uniquify=False)
    if not matches:
        raise chimera.UserError('Could not find any alignment.')
    if not ignore_warnings and len(matches) > 1e6:
        raise chimera.UserError("Too many possible alignments found! This can be "
                                "slow! Use `ignore_warnings true` to try.")
    maps = [list(enumerate(match)) for match in matches]
    best_rmsd, best_xform = 1000., None
    for atom_map in maps:
        rmsd, xform = GetAlignmentTransform(rdk_probe, rdk_reference, atomMap=atom_map,
                                            maxIters=maxIters, reflect=reflect)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_xform = xform

    if transform and best_xform is not None:
        _transform_molecule(probe, chimera_xform(best_xform[:3]))
    return best_rmsd


def align_com(reference, probe, **kwargs):
    ref_com = np.average([a.xformCoord() for a in reference.atoms], axis=0)
    probe_com = np.average([a.xformCoord() for a in probe.atoms], axis=0)
    translation_v = chimera.Point(*ref_com) - chimera.Point(*probe_com)
    xform = probe.openState.xform.translation(translation_v)
    _transform_molecule(probe, xform)
    return 1000.


def cmd_align(reference_sel, probe_sel, methods='best', transform=True, sanitize=True,
              ignore_warnings=False, maxIters=50, reflect=False):
    methods = methods.lower().strip().split(',')
    aligners = []
    for method in methods:
        if method.strip() == 'fast':
            aligners.append(align)
        elif method.strip() == 'fixed':
            aligners.append(untransformed_rmsd)
        elif method.strip() == 'o3a':
            aligners.append(align_o3a)
        elif method.strip() == 'best':
            aligners.append(align_best)
        elif method.strip() == 'com':
            aligners.append(align_com)
        else:
            raise chimera.UserError('{} is not valid for key `method`'.format(method))

    references = reference_sel.molecules()
    probes = probe_sel.molecules()
    if not len(references) == 1:
        raise chimera.UserError("Reference must contain a single molecule.")
    if not len(probes):
        raise chimera.UserError("Select at least one probe.")

    reference = references[0]
    rmsds = []
    for probe in probes:
        if reference.numAtoms < probe.numAtoms:
            raise chimera.UserError("Reference model should be larger than probe.")
        for aligner in aligners:
            try:
                rmsd = aligner(reference, probe, transform=transform, maxIters=maxIters, reflect=reflect,
                                                 sanitize=sanitize, ignore_warnings=ignore_warnings)
            except chimera.UserError as e:
                if aligner is not aligners[-1]:
                    continue
                raise e
        rmsds.append(rmsd)
    msg = ""
    if len(rmsds) == 1:
        msg = "RMSD is {}".format(rmsds[0])
    elif len(rmsds) > 1:
        avg_rmsd = sum(rmsds)/len(rmsds)
        msg = "Average RMSD for {} molecules is {}".format(len(rmsds), avg_rmsd)
    chimera.statusline.show_message(msg, blankAfter=5)
    return rmsds


def cmd_rmsd(reference_sel, probe_sel, method='best', reflect=False, sanitize=False,
             uniquify=False):
    references = reference_sel.molecules()
    probes = probe_sel.molecules()
    if not len(references) == 1:
        raise chimera.UserError("Reference must contain a single molecule.")
    if not len(probes):
        raise chimera.UserError("Select at least one probe.")

    reference = references[0]
    rmsds = []
    for probe in probes:
        if reference.numAtoms < probe.numAtoms:
            raise chimera.UserError("Reference model should be larger than probe.")
        rmsd = untransformed_rmsd(reference, probe, method=method, uniquify=uniquify,
                                  reflect=reflect, sanitize=sanitize)
        rmsds.append(rmsd)
    msg = ""
    if len(rmsds) == 1:
        msg = "RMSD is {}".format(rmsds[0])
    elif len(rmsds) > 1:
        avg_rmsd = sum(rmsds)/len(rmsds)
        print(*rmsds)
        msg = "Average RMSD for {} molecules is {}".format(len(rmsds), avg_rmsd)
    chimera.statusline.show_message(msg, blankAfter=5)
    return rmsds