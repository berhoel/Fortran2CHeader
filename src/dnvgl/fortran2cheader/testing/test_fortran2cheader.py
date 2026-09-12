"""Unit tests for Fortran2CHeader."""

# Copyright (C) 2026 by Berthold Höllmann

import re
from pathlib import Path

import pytest
from dnvgl.fortran2cheader import (
    _ARGS,
    _BIND,
    _FUNCTION,
    _SUBROUTINE,
    _VARTYPE,
    Fortran2CHeader,
)

from .samples import (
    sample_001,
    sample_002,
    sample_003,
    sample_004,
    sample_005,
    sample_006,
)


@pytest.fixture
def h_string_io(tmp_path: Path) -> Path:
    return tmp_path / "test.h"


@pytest.fixture
def pxd_string_io(tmp_path: Path) -> Path:
    return tmp_path / "test.pxd"


SIMPLE_DATA: list[tuple[str, dict[str, str | None]]] = [
    (
        "INTEGER(C_INT), INTENT(IN), VALUE :: iUnit",
        {
            "kind": "C_INT",
            "ftype": "INTEGER",
            "args": "iUnit",
            "modifier": ", INTENT(IN), VALUE ",
            "length": None,
        },
    ),
    (
        "character(kind=c_char), intent(in), dimension(*) :: s",
        {
            "ftype": "character",
            "kind": "c_char",
            "length": None,
            "modifier": ", intent(in), dimension(*) ",
            "args": "s",
        },
    ),
    (
        "character(kind=c_char,len=1), intent(in), dimension(*) :: s",
        {
            "ftype": "character",
            "kind": "c_char",
            "length": "1",
            "modifier": ", intent(in), dimension(*) ",
            "args": "s",
        },
    ),
    (
        "character(kind=c_char,len=1), dimension(*), intent(in) :: s",
        {
            "kind": "c_char",
            "ftype": "character",
            "args": "s",
            "modifier": ", dimension(*), intent(in) ",
            "length": "1",
        },
    ),
    (
        "REAL(C_DOUBLE), INTENT(IN), DIMENSION(n + n + m) :: temp",
        {
            "ftype": "REAL",
            "kind": "C_DOUBLE",
            "length": None,
            "modifier": ", INTENT(IN), DIMENSION(n + n + m) ",
            "args": "temp",
        },
    ),
    (
        "REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: temp",
        {
            "ftype": "REAL",
            "kind": "C_DOUBLE",
            "length": None,
            "modifier": ", INTENT(IN), DIMENSION(n) ",
            "args": "temp",
        },
    ),
]


@pytest.mark.parametrize(argnames=("sample", "reference"), argvalues=SIMPLE_DATA)
def test_simple_1(sample: str, reference: dict[str, str | None]):
    res = _VARTYPE.match(sample)
    assert (res is not None) and res.groupdict() == reference


@pytest.mark.parametrize(
    argnames=("specimen", "expected"), argvalues=[("(s)", {"args": "s"})]
)
def test_args_re_1(specimen, expected):
    probe = re.compile(_ARGS, re.VERBOSE | re.IGNORECASE).match(specimen)
    assert probe is not None
    assert probe.groupdict() == expected


@pytest.mark.parametrize(
    argnames=("specimen", "expected"),
    argvalues=[
        ("bind(c,name='pstr')", {"C": "c", "c_name": "pstr", "quot": "'"}),
        ("BIND(c, name='pstr')", {"C": "c", "c_name": "pstr", "quot": "'"}),
    ],
)
def test_bind_re_1(specimen, expected):
    probe = re.compile(_BIND, re.VERBOSE | re.IGNORECASE).match(specimen)
    assert probe is not None
    assert probe.groupdict() == expected


@pytest.mark.parametrize(
    argnames=("specimen", "expected"),
    argvalues=[
        (
            "subroutine pstr(s) bind(c,name='pstr')",
            {
                "f_name": "pstr",
                "args": "s",
                "C": "c",
                "quot": "'",
                "c_name": "pstr",
            },
        )
    ],
)
def test_subr_re_1(specimen, expected):
    probe = _SUBROUTINE.match(specimen)
    assert probe is not None
    assert probe.groupdict() == expected


@pytest.mark.parametrize(
    argnames=("specimen", "expected"),
    argvalues=[
        (
            "FUNCTION curv2(t, n, x, y, yp, sigma) RESULT(res) BIND(C, NAME='c_curv2')",
            {
                "C": "C",
                "args": "t, n, x, y, yp, sigma",
                "c_name": "c_curv2",
                "f_name": "curv2",
                "prefix": None,
                "quot": "'",
                "result": "res",
            },
        ),
        (
            "FUNCTION pstr(s) RESULT(x) BIND(c, name='pstr')",
            {
                "C": "c",
                "args": "s",
                "c_name": "pstr",
                "f_name": "pstr",
                "prefix": None,
                "quot": "'",
                "result": "x",
            },
        ),
    ],
)
def test_function_re_1(specimen, expected):
    probe = _FUNCTION.match(specimen)
    assert probe is not None
    assert probe.groupdict() == expected


@pytest.mark.parametrize(
    ("i_data", "exp", "exp_pxd"),
    [
        (sample_001.I_DATA, sample_001.EXP, sample_001.EXP_PXD),
        (sample_002.I_DATA, sample_002.EXP, sample_002.EXP_PXD),
        (sample_003.I_DATA, sample_003.EXP, sample_003.EXP_PXD),
        (sample_004.I_DATA, sample_004.EXP, sample_004.EXP_PXD),
        (sample_005.I_DATA, sample_005.EXP, sample_005.EXP_PXD),
        (sample_006.I_DATA, sample_006.EXP, sample_006.EXP_PXD),
    ],
)
def test_subr_2(
    h_string_io: Path, pxd_string_io: Path, i_data: Path, exp: Path, exp_pxd: Path
):
    data = Fortran2CHeader(i_data, signed_to_unsigned_char=True)
    data.parse()
    data.gen_chead(h_string_io)
    res = h_string_io.read_text()
    for i, j in zip(res.split("\n"), exp.read_text().split("\n")):
        if i.startswith("  Generated"):
            continue
        assert i == j
    data.gen_pxd(pxd_string_io)
    res = pxd_string_io.read_text()
    for i, j in zip(res.split("\n"), exp_pxd.read_text().split("\n")):
        if i.startswith("# from "):
            continue
        if i.startswith("# Generated by"):
            continue
        assert i == j
