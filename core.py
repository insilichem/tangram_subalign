#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import cStringIO as StringIO
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform

try:
    from rdkit.Chem import MolFromPDBBlock, SanitizeMol, SanitizeFlags, FastFindRings
    from rdkit.Chem.rdMolAlign import GetAlignmentTransform, GetO3A
    from rdkit.Chem.AllChem import MMFFGetMoleculeProperties, AddHs, GetBestRMS
except ImportError:
    raise chimera.UserError("RDKit must be installed to use this extension.")


def _chimera_to_rdkit(molecule, ignore_metals=True, sanitize=True):
    io = StringIO.StringIO()
    xform = molecule.openState.xform
    if ignore_metals:
        non_metals = [a for a in molecule.atoms if not a.element.isMetal]
        if not non_metals:
            raise ValueError("Molecule does not contain non-metallic atoms.")
        molecule = molecule_from_atoms(molecule, non_metals)
    chimera.pdbWrite([molecule], xform, io)
    io.seek(0)
    rdkit_mol = MolFromPDBBlock(io.getvalue(), False, sanitize)
    io.close()
    if ignore_metals:
        molecule.destroy()
    return rdkit_mol

def _transform_molecule(molecule, xform):
    new_xform = molecule.openState.xform
    new_xform.premultiply(xform)
    molecule.openState.xform = new_xform

def align(reference, probe, transform=True, sanitize=True, **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    rmsd, xform = GetAlignmentTransform(rdk_probe, rdk_reference)
    if transform:
        _transform_molecule(probe, chimera_xform(xform[:3]))
    return rmsd

def align_o3a(reference, probe, transform=True, sanitize=True, **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    FastFindRings(rdk_reference)
    FastFindRings(rdk_probe)
    reference_params = MMFFGetMoleculeProperties(rdk_reference)
    probe_params = MMFFGetMoleculeProperties(rdk_probe)
    o3a = GetO3A(rdk_probe, rdk_reference, probe_params, reference_params)
    rmsd, xform = o3a.Trans()
    if transform:
        _transform_molecule(probe, chimera_xform(xform[:3]))
    return rmsd

def align_best(reference, probe, transform=True, sanitize=True, ignore_warnings=False,
               **kwargs):
    rdk_reference = _chimera_to_rdkit(reference, sanitize=sanitize)
    rdk_probe = _chimera_to_rdkit(probe, sanitize=sanitize)
    matches = rdk_reference.GetSubstructMatches(rdk_probe, uniquify=False)
    if not matches:
        raise chimera.UserError('Could not find any alignment.')
    if ignore_warnings and len(matches) > 1e6:
        raise chimera.UserError("Too many possible alignments found! This can be "
                                "slow! Use `ignore_warnings true` to try.")
    maps = [list(enumerate(match)) for match in matches]
    best_rmsd, best_xform = 1000., None
    for atom_map in maps:
        rmsd, xform = GetAlignmentTransform(rdk_probe, rdk_reference, atomMap=atom_map)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_xform = xform

    if transform and best_xform is not None:
        _transform_molecule(probe, chimera_xform(best_xform[:3]))
    return best_rmsd

def cmd_align(reference_sel, probe_sel, method='best', transform=True, sanitize=True,
              ignore_warnings=False):
    if method.lower().strip() == 'fast':
        aligner = align
    elif method.lower().strip() == 'o3a':
        aligner = align_o3a
    elif method.lower().strip() == 'best':
        aligner = align_best
    else:
        raise ValueError('{} is not valid for key `method`'.format(method))

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
        rmsd = aligner(reference, probe, transform=transform, sanitize=sanitize,
                       ignore_warnings=ignore_warnings)
        rmsds.append(rmsd)
    msg = ""
    if len(rmsds) == 1:
        msg = "RMSD is {}".format(rmsds[0])
    elif len(rmsds) > 1:
        avg_rmsd = sum(rmsds)/len(rmsds)
        msg = "Average RMSD for {} molecules is {}".format(len(rmsds), avg_rmsd)
    chimera.statusline.show_message(msg, blankAfter=5)
