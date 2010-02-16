# $Id: Makefile.am 2248 2006-02-10 15:30:52Z joe $

quadraturerulesdir = $(includedir)/dune/grid/common/quadraturerules
quadraturerules_HEADERS =			\
	cube_imp.hh				\
	gaussquadrature.hh			\
	genericquadrature.hh			\
	jacobi_1_0_imp.hh			\
	jacobi_2_0_imp.hh			\
	lobattoquadrature.hh

noinst_LTLIBRARIES = libquadraturerules.la
libquadraturerules_la_SOURCES = quadraturerules.cc cube.cc jacobi_1_0.cc jacobi_2_0.cc
libquadraturerules_la_LIBADD =			\
	$(DUNE_LIBS)

extra_DIST=jacobian.mw

include $(top_srcdir)/am/global-rules
