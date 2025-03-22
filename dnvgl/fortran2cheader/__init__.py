"""Generate a C/C++ header file from a Fortran source file.

Generate a C/C++ header file from a Fortran source file using those
subroutines decorated with Fortran2003 `C_ISO_BINDINGS` `BIND` and
`C_*` types. `INTERFACE` blocks are ignored to allow direct usage of C
function calls in the `FORTRAN` code.

Parses FORTRAN 90 syntax.

Translation table taken from
<http://de.wikibooks.org/wiki/Fortran:_Fortran_und_C:_Fortran_2003>

+-----------+---------------------------+----------------------+
| Fortran-  | Benannte `iso_c_binding`- | C-Datentyp           |
| Datentyp  | Konstante (kind-Wert)     |                      |
+-----------+---------------------------+----------------------+
|           | c_int                     | int                  |
|           +---------------------------+----------------------+
|           | c_short                   | short int            |
|           +---------------------------+----------------------+
|           | c_long                    | long int             |
|           +---------------------------+----------------------+
|           | c_long_long               | long long int        |
|           +---------------------------+----------------------+
|           | c_signed_char             | signed char,         |
|           |                           | unsigned char        |
|           +---------------------------+----------------------+
|           | c_size_t                  | size_t               |
|           +---------------------------+----------------------+
|           | c_int8_t                  | int8_t               |
|           +---------------------------+----------------------+
|           | c_int16_t                 | int16_t              |
|           +---------------------------+----------------------+
|           | c_int32_t                 | int32_t              |
|           +---------------------------+----------------------+
| integer   | c_int64_t                 | int64_t              |
|           +---------------------------+----------------------+
|           | c_int_least8_t            | int_least8_t         |
|           +---------------------------+----------------------+
|           | c_int_least16_t           | int_least16_t        |
|           +---------------------------+----------------------+
|           | c_int_least32_t           | int_least32_t        |
|           +---------------------------+----------------------+
|           | c_int_least64_t           | int_least64_t        |
|           +---------------------------+----------------------+
|           | c_int_fast8_t             | int_fast8_t          |
|           +---------------------------+----------------------+
|           | c_int_fast16_t            | int_fast16_t         |
|           +---------------------------+----------------------+
|           | c_int_fast32_t            | int_fast32_t         |
|           +---------------------------+----------------------+
|           | c_int_fast64_t            | int_fast64_t         |
|           +---------------------------+----------------------+
|           | c_intmax_t                | intmax_t             |
|           +---------------------------+----------------------+
|           | c_intptr_t                | intptr_t             |
+-----------+---------------------------+----------------------+
|           | c_float                   | float                |
|           +---------------------------+----------------------+
| real      | c_double                  | double               |
|           +---------------------------+----------------------+
|           | c_long_double             | long double          |
+-----------+---------------------------+----------------------+
|           | c_float_complex           | float _Complex       |
|           +---------------------------+----------------------+
| complex   | c_double_complex          | double _Complex      |
|           +---------------------------+----------------------+
|           | c_long_double_complex     | long double _Complex |
+-----------+---------------------------+----------------------+
| logical   | c_bool                    | _Bool                |
+-----------+---------------------------+----------------------+
| character | c_char                    | char                 |
+-----------+---------------------------+----------------------+
"""

from __future__ import annotations

import argparse
import os
import os.path
from pathlib import Path
import re
import sys
import types
from typing import TYPE_CHECKING, ClassVar

from rich.console import Console

if TYPE_CHECKING:
    from collections.abc import Generator

__date__ = "0000/00/00 00:00:00 xxx"
__author__ = "`Berthold Höllmann <berthold.hoellmann@dnvgl.com>`__"
__copyright__ = "Copyright © 2010 by DNV GL SE"
__credits__ = ["Berthold Höllmann"]
__maintainer__ = "Berthold Höllmann"
__email__ = "berhoel@gmail.com"

_CONSOLE = Console()


def casi(inp: str) -> str:
    """Make string for case insenitive regular expression string.

    Example::

      >>> casi('abc') == '[aA][bB][cC]'
      True
    """
    outp: list[str] = []
    letter = re.compile(r"^\w$")
    whitespace = re.compile(r"^\s$")
    for char in inp:
        if letter.match(char):
            outp.extend(("[", char.lower(), char.upper(), "]"))
        elif whitespace.match(char):
            outp.append("\\s")
        else:
            outp.append(char)
    return "".join(outp)


# SUBROUTINE SXFGeRh (iUnit, oRelName, oNoOAttr, oNoORows, oAttName,
#                     oAttType, oAttLeng) BIND(C,NAME="SXFGeRh")
_BIND = rf"""{casi("BIND")}
\( \s*
(?:
  (?:
    (?P<C>[Cc]) |
    (?:{casi("NAME")} ) \s* =
       \s* (?P<quot>[\'\"]) (?P<c_name>[\w]+) (?P=quot)
  ) \s* ,? \s*
)+
\s* \)
"""
# The first solution becomes to slow to long argument lists.
_ARGS = (
    r"\( \s* (?P<args> (?: \w+? \s* ,? \s* )* ) \s* \) \s*"
    if False
    else r"\( \s* (?P<args> .*? ) \s* \) \s*"
)
_SUBROUTINE = re.compile(
    rf"""
^
(?: {casi("SUBROUTINE")} ) \s+
(?P<f_name> \w+ ) \s*
{_ARGS}{_BIND}""",
    re.VERBOSE,
)

# FUNCTION C_CALLOC(elt_count, elt_size) RESULT(ptr) BIND(C, NAME="calloc")
_FUNCTION = re.compile(
    r"""
^
(?P<prefix> .+)?? \s*
(?: {casi("FUNCTION")} ) \s+
(?P<f_name> \w+ ) \s*
{_ARGS}
(?: {casi("RESULT")} \s*
[(] \s* (?P<result> .+? ) [)] ) \s*
{_BIND}""",
    re.VERBOSE,
)
# INTEGER(C_INT), INTENT(IN), VALUE :: iUnit
_VARTYPE = re.compile(
    rf"""
^
(?P<ftype>
  (?: {casi("INTEGER")} ) |
  (?: {casi("REAL")} ) |
  (?: {casi("COMPLEX")} ) |
  (?: {casi("LOGICAL")} ) |
  (?: {casi("CHARACTER")} ) |
  (?: {casi("TYPE")} )
)
\( \s*
  (?:
    (?:
      (?: (?: {casi("KIND=")} )? (?P<kind> [\w\d]+ ) )? |
      (?: (?: {casi("LEN=")} ) (?P<length> [\d]+ ) )?
    ) \s* ,?
  )*
\s* \) \s*
(?P<modifier> (?: , \s* [*\w()]+ \s* )+ )? :: \s*
(?P<args> (?: [\w]+ \s* ,? \s* )* )
""",
    re.VERBOSE,
)

_INTENT = re.compile(
    rf""".*{casi("intent")}\s*\(\s*(?P<dir>
(?:{casi("IN")})|
(?:{casi("OUT")})|
(?:{casi("INOUT")})|
(?:{casi("IN,OUT")}))\s*\).*""",
    re.VERBOSE,
)

_TYPE = re.compile(
    rf"""(?P<ftype>
  (?: {casi("INTEGER")} ) |
  (?: {casi("REAL")} ) |
  (?: {casi("COMPLEX")} ) |
  (?: {casi("LOGICAL")} ) |
  (?: {casi("CHARACTER")} ) |
  (?: {casi("TYPE")} )
)
\( \s* (?P<kind> [\w\d=]+ ) \s* \)""",
    re.VERBOSE,
)

_INTERFACE = re.compile(rf"^{casi('INTERFACE')}$")
_END_INTERFACE = re.compile(rf"^{casi('END INTERFACE')}$")


def file_newer(new: Path, old: Path) -> bool:
    """Return if `new` is newer than `old`."""
    n_time = new.stat().st_mtime
    o_time = old.stat().st_mtime
    return n_time > o_time


class FortranSourceProvider:
    """Provide concatenated Fortran source lines for analysis."""

    def __init__(self, data: Path) -> None:
        """Initialize class instance."""
        self.file = data
        self._read = self.__read()

    def __iter__(self) -> FortranSourceProvider:
        """Return iterable."""
        return self

    def __read(self) -> Generator:
        """Return next raw line from input file.

        Ignoring empty lines and comments.
        """
        with self.file.open() as inp:
            for line in inp:
                pos = line.find("!")
                if pos >= 0:
                    yield line[:pos].strip()
                else:
                    yield line.strip()

    def next(self) -> str:
        """Return next line."""
        line = next(self._read)
        while line.endswith("&"):
            line = line[:-1]
            n_line = next(self._read)
            n_line = n_line.removeprefix("&")
            line += n_line
        return line.strip()

    __next__ = next


class Comment:
    """Provide C comments."""

    flavour = "C"

    def __init__(self, text: str) -> None:
        """Initialize class instance."""
        self.text = text

    def __str__(self) -> str:
        """Return string representation."""
        text = self.text.split("\n")
        if Comment.flavour == "C":
            content = "\n".join((f"  {i}").rstrip() for i in text)
            return f"""/*
{content}
 */
"""
        if Comment.flavour == "pxd":
            return "\n".join((f"# {t}").strip() for t in text)
        return "\n".join(text)


class _Routine:
    """Base class for representing Fortran routines."""

    f_kinds: ClassVar[dict[str, dict[str, str]]] = {
        "integer": {
            "c_int": "int",
            "c_short": "short int",
            "c_long": "long int",
            "c_long_long": "long long int",
            "c_size_t": "size_t",
            "c_int8_t": "int8_t",
            "c_int16_t": "int16_t",
            "c_int32_t": "int32_t",
            "c_int64_t": "int64_t",
            "c_int_least8_t": "int_least8_t",
            "c_int_least16_t": "int_least16_t",
            "c_int_least32_t": "int_least32_t",
            "c_int_least64_t": "int_least64_t",
            "c_int_fast8_t": "int_fast8_t",
            "c_int_fast16_t": "int_fast16_t",
            "c_int_fast32_t": "int_fast32_t",
            "c_int_fast64_t": "int_fast64_t",
            "c_intmax_t": "intmax_t",
            "c_intptr_t": "intptr_t",
        },
        "real": {
            "c_float": "float",
            "c_double": "double",
            "c_long_double": "long double",
        },
        "complex": {
            "c_float_complex": "float _Complex",
            "c_double_complex": "double _Complex",
            "c_long_double_complex": "long double _Complex",
        },
        "logical": {"c_bool": "char"},  # "_Bool"},
        "character": {"c_char": "char"},
        "type": {"c_ptr": "void*", "c_funptr": "(*)"},
    }

    def __init__(self, *, signed_to_unsigned_char: bool):
        """Intitialize class instance."""
        if signed_to_unsigned_char:
            self.f_kinds["integer"]["c_signed_char"] = "unsigned char"
        else:
            self.f_kinds["integer"]["c_signed_char"] = "signed char"
        self.argdict: dict[str, list[None | str | Comment]] = {}
        self.uargs: list[str] = []
        self.comment: Comment | None = None
        self.name: str | None = None
        self.result: str | None = None

    def __str__(self) -> str:
        """Return string representation."""
        out = f"{self.comment}"
        args = ", ".join(
            f"{self.argdict[arg][0]} {self.argdict[arg][1]}" for arg in self.uargs
        )
        if Comment.flavour == "C":
            if not args:
                args = "void"
            proto = f"extern {self.result} {self.name}({args});\n"
        elif Comment.flavour == "pxd":
            proto = f"""
    {self.result} {self.name}({args.replace("_Complex", "complex")})"""

        out += proto
        return out

    def add_arg(
        self, args: str, ftype: str, kind: str, modifier: str, length: str
    ) -> str | None:
        """Add argument information to Subroutine information."""
        _ = length
        c_type = self.f_kinds.get(ftype.lower(), {}).get(kind.lower(), None)
        intent = modifier and _INTENT.match(modifier)
        if c_type and modifier and intent and intent.groupdict()["dir"].lower() == "in":
            c_type = "const " + c_type
        if (
            c_type
            and modifier
            and ("value" not in modifier.lower() or "dimension" in modifier.lower())
        ):
            c_type += "*"
        for arg in (a.strip().upper() for a in args.split(",")):
            if arg in self.uargs:
                self.argdict[arg][0] = c_type
        return c_type


class Subroutine(_Routine):
    """Representing Fortran SUBBROUTINEs."""

    def __init__(
        self,
        line: str,
        c_name: str,
        f_name: str,
        args: list[str],
        *,
        signed_to_unsigned_char: bool,
    ):
        """Initialize class instance."""
        super().__init__(signed_to_unsigned_char=signed_to_unsigned_char)
        self.comment = Comment(
            f"""{c_name}
Generated from FORTRAN routine '{f_name}'
FORTRAN declaration:
    {line}"""
        )
        self.name = c_name
        args = [a.strip() for a in args if a]
        self.uargs = [a.upper() for a in args]
        self.argdict = {a.upper(): [None, a] for a in args}
        self.result = "void"


class Function(_Routine):
    """Representing Fortran FUNCTIONs."""

    def __init__(  # noqa:PLR0913
        self,
        line: str,
        c_name: str,
        f_name: str,
        prefix: str,
        result: str,
        args: list[str],
        *,
        signed_to_unsigned_char: bool,
    ) -> None:
        """Initialize class instance."""
        super().__init__(signed_to_unsigned_char=signed_to_unsigned_char)
        self.comment = Comment(
            f"""\
{c_name}
Generated from FORTRAN routine '{f_name}'
FORTRAN declaration:
    {line}"""
        )

        self.name = c_name
        args = [a.strip() for a in args if a]
        self.uargs = [a.upper() for a in args]
        self.argdict = {a.upper(): [None, a] for a in args}
        if result:
            self.resultN = result
        else:
            self.resultN = f_name
        if _prefix := _TYPE.match(prefix):
            self.result = self.f_kinds.get(_prefix.group("ftype").lower(), {}).get(
                _prefix.group("kind").lower(), None
            )
        else:
            self.result = None

    def add_arg(
        self, args: str, ftype: str, kind: str, modifier: str, length: str
    ) -> None:
        """Add argument information to Subroutine information."""
        c_type = super().add_arg(args, ftype, kind, modifier, length)
        if self.resultN.upper() in [a.strip().upper() for a in args.split(",")]:
            self.result = c_type


class Fortran2CHeader:
    """Extract a C header file from a Fortran file.

    Extract using `BIND(C,name='xxx')` for providing a C compatible interface.

    Only arguments with type kinds from `ISO_C_BINDING` module.
    """

    def __init__(self, data: Path, **kw: bool) -> None:
        """Initialize Class."""
        self.input = data
        self.name = data.name
        self.basename = data.stem

        self.data = FortranSourceProvider(self.input)

        self.signed_to_unsigned_char = kw.get("signed_to_unsigned_char", False)
        self.force = kw.get("force", False)
        self.generate_pxd = kw.get("generate_pxd", True)

    def parse(self) -> None:  # noqa:PLR0912,C901
        """Parse the input file for `ISO_C_BINDING` information."""
        fname = ""

        if isinstance(self.input, types.GeneratorType):
            fname = "generator"
        else:
            fname = self.input.name
        _CONSOLE.print(f"*** fortran2cheader - Parsing {fname}")
        subr: None | _Routine = None
        interface = False
        self.info = []
        for i in self.data:
            vartype = _VARTYPE.match(i)
            if not interface and (line := _SUBROUTINE.match(i)):
                if subr:
                    self.info.append(subr)
                gdict = line.groupdict()
                if not gdict["C"]:
                    subr = None
                    continue
                subr = Subroutine(
                    signed_to_unsigned_char=self.signed_to_unsigned_char,
                    f_name=gdict["f_name"],
                    c_name=gdict["c_name"],
                    args=gdict["args"].split(","),
                    line=i,
                )
            elif not interface and (line := _FUNCTION.match(i)):
                if subr:
                    self.info.append(subr)
                gdict = line.groupdict()
                if not gdict["C"]:
                    subr = None
                    continue
                subr = Function(
                    signed_to_unsigned_char=self.signed_to_unsigned_char,
                    f_name=gdict["f_name"],
                    c_name=gdict["c_name"],
                    args=gdict["args"].split(","),
                    prefix=gdict["prefix"],
                    result=gdict["result"],
                    line=i,
                )
            elif _INTERFACE.match(i):
                interface = True
            elif _END_INTERFACE.match(i):
                interface = False
            elif not interface and subr and vartype and vartype.groupdict()["kind"]:
                gdict = vartype.groupdict()
                subr.add_arg(
                    args=gdict["args"],
                    ftype=gdict["ftype"],
                    kind=gdict["kind"],
                    modifier=gdict["modifier"],
                    length=gdict["length"],
                )

        if subr:
            self.info.append(subr)

    def gen_chead(self, outf: Path) -> None:
        """Generate the output file."""
        Comment.flavour = "C"
        with outf.open("w") as sink:
            sink.write(
                "\n".join(
                    f"{s}"
                    for s in (
                        Comment(
                            f"""\
{outf.name}
Header file generated from parsing ISO_C_BINDING information
from {self.name}.

Generated by {sys.argv[0]}, version {
                                __import__(
                                    "importlib.metadata", fromlist=["version"]
                                ).version("Fortran2CHeader")
                            }."""
                        ),
                        f"""\
#ifndef {self.basename.upper()}_H
#define {self.basename.upper()}_H

#ifdef __cplusplus
extern "C" {{
#endif /* __cplusplus */

""",
                    )
                )
            )
            sink.write("\n".join(f"{s}" for s in self.info))
            sink.write(
                f"""
#ifdef __cplusplus
}} /* extern "C" */
#endif /* __cplusplus */

#endif /* {self.basename.upper()}_H */"""
            )

    def gen_pxd(self, outf: Path, header: str | None = None) -> None:
        """Generate the output file."""
        if header is None:
            header = f"{self.basename}.h"
        Comment.flavour = "pxd"
        with outf.open("w") as sink:
            sink.write(
                "\n".join(
                    f"{s}".rstrip()
                    for s in (
                        Comment(
                            f"""
{outf.name}
Cython Header file generated from parsing ISO_C_BINDING information
from {self.input}.

Generated by {os.path.split(sys.argv[0])[-1]}, version {
                                __import__(
                                    "importlib.metadata", fromlist=["version"]
                                ).version("Fortran2CHeader")
                            }.""",
                        ),
                        "",
                        f'cdef extern from "{header}" nogil:',
                        "",
                    )
                )
            )
            sink.write("\n".join((f"{s}").rstrip() for s in self.info))

    def gen_output(self, h_name: Path, pxd_name: Path) -> None:
        """Generate output."""
        if not self.force and h_name.is_file() and file_newer(h_name, self.input):
            _CONSOLE.print(
                "*** fortran2cheader - output file is newer than input, "
                f"keeping '{h_name}'."
            )
        else:
            _CONSOLE.print(f"*** fortran2cheader - generating output '{h_name}'.")
            self.gen_chead(Path(h_name))

        if not self.generate_pxd:
            return

        if not self.force and pxd_name.exists() and file_newer(pxd_name, self.input):
            _CONSOLE.print(f"*** fortran2cheader - generating output '{pxd_name}'.")
        else:
            _CONSOLE.print(
                "*** fortran2cheader - output file is newer than input, "
                f"keeping '{pxd_name}'."
            )
            self.gen_pxd(Path(pxd_name))


def parse_cmdline() -> argparse.ArgumentParser:
    """Parse the command line options."""
    parser = argparse.ArgumentParser(
        description="""
Generate a C/C++ header file from a Fortran file using C_ISO_BINDINGS."""
    )
    parser.add_argument("infile", type=Path, help="Fortran input file")
    parser.add_argument(
        "--signed-to-unsigned-char",
        "-s",
        action="store_true",
        default=False,
        help="""
Use 'unsigned char' instead for 'signed char' for 'c_signed_char'""",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        default=False,
        help="""
Force generation of output files even if they exist, and are newer
than the input file.""",
    )
    parser.add_argument(
        "--generate-pxd",
        "-p",
        action="store_true",
        default=False,
        help="""
Generate also a pxd file for import in Cython process.""",
    )
    return parser


class Fortran2CHeaderCMD(Fortran2CHeader):
    """Command line interface for Fortran2CHeader."""

    def __init__(self) -> None:
        """Initialize Class."""
        options = parse_cmdline().parse_args()
        super().__init__(
            data=options.infile,
            signed_to_unsigned_char=options.signed_to_unsigned_char,
            force=options.force,
            generate_pxd=options.generate_pxd,
        )


def main() -> None:
    """Execute main program."""
    f_file = Fortran2CHeaderCMD()

    h_name = Path(f"{f_file.basename}.h")
    pxd_name = Path(f"{f_file.basename}.pxd")

    f_file.parse()
    f_file.gen_output(h_name, pxd_name)


if __name__ == "__main__":
    main()
