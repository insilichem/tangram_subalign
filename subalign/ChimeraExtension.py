#!/usr/bin/env python
# -*- coding: utf-8 -*-


# import chimera.extension
from Midas.midas_text import addCommand

# class SubAlignEMO(chimera.extension.EMO):

#     def name(self):
#         return 'SubAlign'

#     def description(self):
#         return 'Substructure Alignment'

#     def categories(self):
#         return ['InsiliChem']

#     def icon(self):
#         pass

#     def activate(self):
#         pass

# chimera.extension.manager.registerExtension(SubAlignEMO(__file__))


def cmd_subalign(cmdName, args):
    from Midas.midas_text import doExtensionFunc
    from subalign.core import cmd_align
    doExtensionFunc(cmd_align, args, specInfo=[("refSpec", "reference_sel", None),
                                               ("probeSpec", "probe_sel", None)])

def cmd_subrmsd(cmdName, args):
    from Midas.midas_text import doExtensionFunc
    from subalign.core import cmd_rmsd
    doExtensionFunc(cmd_rmsd, args, specInfo=[("refSpec", "reference_sel", None),
                                               ("probeSpec", "probe_sel", None)])

addCommand("subalign", cmd_subalign)
addCommand("subrmsd", cmd_subrmsd)
