# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= poetry run sphinx-build
SPHINXAPIDOC  ?= poetry run sphinx-apidoc
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile API-GEN
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

API-GEN:
	$(SPHINXAPIDOC)			\
	 --module-first			\
	 --force				\
	 --private				\
	 --implicit-namespaces			\
	 --doc-project Fortran2CHeader          \
	 -o _tmp				\
	 ../dnvgl				\
	 ../dnvgl/fortran2cheader/testing/

publish: html
	lftp software_docu -e "mirror -R _build/html/ python/Fortran2CHeader;exit"
	curl "$(shell idn https://www.google.com/ping?sitemap=https://python.höllmanns.de/Fortran2CHeader/sitemap.xml)"

