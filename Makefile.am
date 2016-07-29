# Copyright 2016, Alistair Boyle, 3-clause BSD License
# build to GNU standards
AUTOMAKE_OPTIONS = foreign

# the actual program to be installed at the end
bin_PROGRAMS = zedhat
zedhat_SOURCES = src/zedhat.c src/argv.c

#AM_CFLAGS = $(CFLAGS) -Werror -Wall #-Wextra
AM_CFLAGS = $(CFLAGS) -Werror -Wall #-Wextra
LIBS += $(LAPACKE_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# header files that need to be distributed but not installed
noinst_HEADERS = src/argv.h

# want these scripts installed in the end
#bin_SCRIPTS = 

# some random files we want to distribute with the package tar
EXTRA_DIST = $(bin_SCRIPTS) \
             ChangeLog README.md NEWS LICENSE

BUILT_TESTS = tests/libcblas_dgemm tests/liblapacke_dgels tests/libumfpack tests/libmatio
check_PROGRAMS = $(BUILT_TESTS)
#tests_unit_perftimer_SOURCES = tests/unit-perftimer.c src/perftimer.c
#tests_unit_perftimer_CPPFLAGS = -I$(srcdir)/src

clean-local: clean-local-check
	  -rm -f src/*.lo

maintainer-clean-local:
	-rm -f aclocal.m4 configure ChangeLog Makefile.in src/config.h.in*
	-rm -rf build-aux .deps

# regenerate the ChangeLog from git's history via ./autogen (updates configure, etc too)
# don't fail if autogen doesn't exist though: its not in the .tar-ed releases
dist-hook: $(top_srcdir)/autogen
	[[ -x $(top_srcdir)/autogen ]] && $(top_srcdir)/autogen $(distdir) $(top_srcdir) || true

.PHONY: cov
cov: dist
	$(srcdir)/cov $(distdir).tar.gz

# apply astyle to source files
style:
	astyle --suffix=none --options=./astylerc src/*.{c,h} tests/*.c

commit:
	make clean && CC=clang ./configure && make check
	make clean && ./configure && make check
	make style
	git ls-files -m
	[ "x$(shell git ls-files -m | wc -l)" == "x0" ] # no uncommitted files?

######
## vodoo incantations to get the test suite automation to build (autotest)

# The `:;' works around a Bash 3.2 bug when the output is not writeable.
$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	:;{ \
	  echo '# Signature of the current package.' && \
	  echo 'm4_define([AT_PACKAGE_NAME],' && \
	  echo '  [$(PACKAGE_NAME)])' && \
	  echo 'm4_define([AT_PACKAGE_TARNAME],' && \
	  echo '  [$(PACKAGE_TARNAME)])' && \
	  echo 'm4_define([AT_PACKAGE_VERSION],' && \
	  echo '  [$(PACKAGE_VERSION)])' && \
	  echo 'm4_define([AT_PACKAGE_STRING],' && \
	  echo '  [$(PACKAGE_STRING)])' && \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT],' && \
	  echo '  [$(PACKAGE_BUGREPORT)])'; \
	  echo 'm4_define([AT_PACKAGE_URL],' && \
	  echo '  [$(PACKAGE_URL)])'; \
	} >'$(srcdir)/package.m4'


EXTRA_DIST += testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in
TESTSUITE = $(srcdir)/testsuite
test: check

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
	  $(TESTSUITEFLAGS) -k '!hello,!unit'
# don't run dependancy and unit tests

clean-local-check:
	test ! -f '$(TESTSUITE)' || \
	  $(SHELL) '$(TESTSUITE)' --clean
	-rm -f $(TESTSUITE)
	-rm -f $(srcdir)/package.m4
	-rm -f atconfig

#atconfig: $(srcdir)/config.status
#	$(srcdir)/config.status

AUTOM4TE = $(SHELL) $(abs_srcdir)/build-aux/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(abs_srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@


# valgrind configuration
#VALGRIND_BIN=valgrind
#VALGRIND_ENV=MPIWRAP_DEBUG=${MPIWRAP_DEBUG_ARGS} LD_PRELOAD=${MPI_PREFIX}/lib/valgrind/libmpiwrap-${MPI_PLATFORM}.so
#VALGRIND_ARGS=--trace-children=yes --track-origins=yes --leak-check=full
#--db-attach=yes
#MPI_BIN=mpirun
#MPI_ARGS=-n 4
#VALGRIND=${VALGRIND_ENV} ${VALGRIND_BIN} ${VALGRIND_ARGS} ${MPI_BIN} ${MPI_ARGS}
#valgrind: tests/mumps-file
#	cd tests && ${VALGRIND} mumps-file

.PHONY: help
help:
	@echo "Useful make targets: make [target]"
	@echo "  <none> == all"
	@echo "  all              - build software"
	@echo "  check            - run testsuite"
	@echo "  install"
	@echo "  uninstall"
	@echo "  clean            - remove built files"
	@echo "For developers:"
	@echo "  distclean        - clean to distributed state"
	@echo "  maintainer-clean - delete all generated files, requires ./autogen to rebuild"
	@echo "  dist             - create a .tar.gz for distribution"
	@echo "  distcheck        - .tar, then untar in new dir and test"
	@echo "  cov              - report test coverage"
	@echo "  tags             - generate ctags for editing source files"
	@echo "  style            - run astyle on source files to tidy formatting"
