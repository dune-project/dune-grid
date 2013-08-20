CURDIR = doc/latex
BASEDIR = ../..

docdir=$(datadir)/doc/dune-grid

doc_DATA =
EXTRA_DIST =
CLEANFILES =

doc_DATA += hybrid-grid.pdf
EXTRA_DIST += hybrid-grid.tex hybrid-grid.bib
hybrid-grid.dvi: hybrid-grid.tex hybrid-grid.bib
CLEANFILES += hybrid-grid.pdf

EXTRA_DIST += CMakeLists.txt

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/latex
