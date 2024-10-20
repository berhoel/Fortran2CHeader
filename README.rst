Fortran2CHeader
===============

Generate a C/C++ header file from a Fortran source file using those
subroutines decorated with Fortran2003 ``C_ISO_BINDINGS`` ``BIND`` and
``C_*`` types. ``INTERFACE`` blocks are ignored to allow direct usage of
C function calls in the ``FORTRAN`` code.

Python based tool to generate Cython ``pxd`` and C/C++ header file from
Fortran file using ``ISO_C_BINDINGS``. To user from Python ``setup.py``
file use

.. code:: python

   ...
   from pathlib import Path
   from dnvgl.fortran2cheader import Fortran2CHeader

   HEADER = Fortran2CHeader(
       data=open(os.path.join('xx.f90')),
       signed_to_unsigned_char=False)
   HEADER.parse()
   CHEAD_NAME = 'xx.h'
   PXD_NAME = 'xx.pxd'
   HEADER.gen_output(Path(CHEAD_NAME), Path(PXD_NAME))
   ...

Then use ``xx.pxd`` in your Cython module.
