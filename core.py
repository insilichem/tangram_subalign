#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import cStringIO as StringIO
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform

try:
    from rdkit.Chem import MolFromPDBBlock
    from rdkit.Chem.rdMolAlign import GetAlignmentTransform
except ImportError:
    raise chimera.UserError("RDKit must be installed to use this extension.")


def _chimera_to_rdkit(molecule, ignore_metals=True):
    io = StringIO.StringIO()
    xform = molecule.openState.xform
    if ignore_metals:
        non_metals = [a for a in molecule.atoms if not a.element.isMetal]
        if not non_metals:
            raise ValueError("Molecule does not contain non-metallic atoms.")
        molecule = molecule_from_atoms(molecule, non_metals)
    chimera.pdbWrite([molecule], xform, io)
    io.seek(0)
    rdkit_mol = MolFromPDBBlock(io.getvalue(), False, False)
    io.close()
    if ignore_metals:
        molecule.destroy()
    return rdkit_mol


def align(reference, probe, transform=True):
    rdk_reference = _chimera_to_rdkit(reference)
    rdk_probe = _chimera_to_rdkit(probe)
    rmsd, xform = GetAlignmentTransform(rdk_probe, rdk_reference)
    if transform:
        new_xform = probe.openState.xform
        new_xform.premultiply(chimera_xform(xform[:3]))
        probe.openState.xform = new_xform
    return rmsd


def cmd_align(reference_sel, probe_sel, transform=True):
    references = reference_sel.molecules()
    probes = probe_sel.molecules()
    if not len(references) == 1:
        raise chimera.UserError("Reference must contain a single molecule.")
    if not len(probes):
        raise chimera.UserError("Select at least one probe.")

    reference = references[0]
    rmsds = []
    for probe in probes:
        rmsd = align(reference, probe, transform=transform)
        rmsds.append(rmsd)
    msg = ""
    if len(rmsds) == 1:
        msg = "RMSD is {}".format(rmsds[0])
    elif len(rmsds) > 1:
        avg_rmsd = sum(rmsds)/len(rmsds)
        msg = "Average RMSD for {} molecules is {}".format(len(rmsds), avg_rmsd)
    chimera.statusline.show_message(msg, blankAfter=5)
        

