# which tests to run
TESTS =

# programs just to build when "make check" is used
check_PROGRAMS =

# define the programs
TESTS += structuredgridfactorytest
check_PROGRAMS += structuredgridfactorytest
structuredgridfactorytest_SOURCES = structuredgridfactorytest.cc
structuredgridfactorytest_CPPFLAGS = $(AM_CPPFLAGS) \
	                             $(ALUGRID_CPPFLAGS)			\
	                             $(UG_CPPFLAGS)
structuredgridfactorytest_LDFLAGS = $(AM_LDFLAGS) \
	                            $(ALUGRID_CPPFLAGS)			\
	                            $(UG_LDFLAGS)
structuredgridfactorytest_LDADD = $(LDADD) \
	                                $(ALUGRID_LIBS)				\
                                  $(UG_LIBS)

TESTS += vertexordertest
check_PROGRAMS += vertexordertest
vertexordertest_SOURCES = vertexordertest.cc
vertexordertest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)			\
	$(UG_CPPFLAGS)
vertexordertest_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALUGRID_LDFLAGS)			\
	$(UG_LDFLAGS)
vertexordertest_LDADD =				\
	$(UG_LIBS)				\
	$(ALUGRID_LIBS)				\
	$(LDADD)

TESTS += persistentcontainertest
check_PROGRAMS += persistentcontainertest
persistentcontainertest_SOURCES = persistentcontainertest.cc
persistentcontainertest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)			\
	$(UG_CPPFLAGS)
persistentcontainertest_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALUGRID_LDFLAGS)			\
	$(UG_LDFLAGS)
persistentcontainertest_LDADD =				\
	$(UG_LIBS)				\
	$(ALUGRID_LIBS)				\
	$(LDADD)

include $(top_srcdir)/am/global-rules
