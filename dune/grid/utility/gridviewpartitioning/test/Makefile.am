check_PROGRAMS =
check_SCRIPTS =
TESTS =
CLEANFILES =

AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS)
LDADD =						\
	$(DUNEMPILIBS)				\
	@LDADD@
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)

check_HEADERS =					\
	evolve.hh				\
	finitevolume_main.hh			\
	initialize.hh				\
	transportproblem2.hh			\
	vtkout.hh

if ALUGRID

TESTS += finitevolume_alusimplex.sh
check_SCRIPTS += finitevolume_alusimplex.sh
check_PROGRAMS += finitevolume_alusimplex
finitevolume_alusimplex_SOURCES = finitevolume_alusimplex.cc
finitevolume_alusimplex_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
finitevolume_alusimplex_LDADD =			\
	$(ALUGRID_LIBS)				\
	$(LDADD)
finitevolume_alusimplex_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
CLEANFILES +=					\
	alusimplex_*_concentration-*.vtu	\
	alusimplex_*_concentration.series

endif ALUGRID

TESTS += finitevolume_yasp.sh
check_SCRIPTS += finitevolume_yasp.sh
check_PROGRAMS += finitevolume_yasp
finitevolume_yasp_SOURCES = finitevolume_yasp.cc
CLEANFILES +=					\
	yasp_*_concentration-*.vtu		\
	yasp_*_concentration.series

include $(top_srcdir)/am/global-rules
