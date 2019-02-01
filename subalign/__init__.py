#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from .core import align, align_o3a, align_best, untransformed_rmsd
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
