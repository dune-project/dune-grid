# which tests to run
TESTS =

# programs just to build when "make check" is used
check_PROGRAMS =

# define the programs
TESTS += structuredgridfactorytest
check_PROGRAMS += structuredgridfactorytest
structuredgridfactorytest_SOURCES = structuredgridfactorytest.cc
structuredgridfactorytest_CPPFLAGS = $(AM_CPPFLAGS) \
	                             $(UG_CPPFLAGS)
structuredgridfactorytest_LDFLAGS = $(AM_LDFLAGS) \
	                            $(UG_LDFLAGS)
structuredgridfactorytest_LDADD = $(LDADD) \
                                  $(UG_LIBS)

include $(top_srcdir)/am/global-rules
