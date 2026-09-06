"""Provide set of sampple data."""

# Copyright (C) 2026 by Berthold Höllmann

from pathlib import Path

from . import HERE

NAME = Path(__file__).stem

I_DATA = HERE / f"{NAME}.f90"
EXP = HERE / f"{NAME}.h"
EXP_PXD = HERE / f"{NAME}.pxd"
