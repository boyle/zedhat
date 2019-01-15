# Copyright 2016, Alistair Boyle, 3-clause BSD License
# build to GNU standards
AUTOMAKE_OPTIONS = foreign

# the actual program to be installed at the end
bin_PROGRAMS = zedhat
zedhat_SOURCES = src/zedhat.c src/argv.c src/matrix.c

AM_CFLAGS = $(CFLAGS)
AM_CFLAGS += $(MATIO_CFLAGS) $(UMFPACK_CFLAGS) $(LAPACKE_CFLAGS) $(LAPACK_CFLAGS) $(BLAS_CFLAGS)
LIBS += $(MATIO_LIBS) $(UMFPACK_LIBS) $(CHOLMOD_LIBS) $(LAPACKE_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
AM_CFLAGS += -Werror -Wall #-Wextra

# header files that need to be distributed but not installed
noinst_HEADERS = src/argv.h src/matrix.h src/model.h src/readngvol.h

dist_man_MANS = man/zedhat.1

# want these scripts installed in the end
#bin_SCRIPTS =

# some random files we want to distribute with the package tar
EXTRA_DIST = $(bin_SCRIPTS) \
				 ChangeLog README.md NEWS LICENSE

BUILT_TESTS = tests/libcblas_dgemm tests/liblapacke_dgels \
			  tests/libumfpack tests/libcholmod tests/libmatio \
			  tests/readngvol tests/model
check_PROGRAMS = $(BUILT_TESTS)
#tests_unit_perftimer_SOURCES = tests/unit-perftimer.c src/perftimer.c
#tests_unit_perftimer_CPPFLAGS = -I$(srcdir)/src
#readngvol_SOURCES = src/readngvol.c
tests_readngvol_SOURCES = tests/readngvol.c src/readngvol.c src/model.c
tests_model_SOURCES = tests/model.c src/model.c tests/cmocka.c

clean-local: clean-local-check
	  -rm -f src/*.lo

maintainer-clean-local:
	-rm -f aclocal.m4 configure ChangeLog Makefile.in src/config.h.in*
	-rm -rf build-aux .deps

# regenerate the ChangeLog from git's history via ./autogen (updates configure, etc too)
# don't fail if autogen doesn't exist though: its not in the .tar-ed releases
dist-hook:
#	$(top_srcdir)/autogen
	[[ -x $(top_srcdir)/autogen ]] && $(top_srcdir)/autogen $(distdir) $(top_srcdir) || true
	-rm -rf $(distdir)/autom4te.cache

.PHONY: cov
cov.log: dist
	$(srcdir)/cov | tee cov.log
cov: cov.log
#	[ "x$(shell tail -4 $^ | grep '100.0%' | wc -l)" == "x3" ] || exit 1 # 100% test coverage please!

STYLE_FILES=src/*.c src/*.h tests/*.c
style: style.log
style.log: $(STYLE_FILES) astylerc
	astyle --suffix=none --options=./astylerc $(STYLE_FILES) | tee $@
check_style: style.log
	-grep Formatted $^
	[ "x$(shell grep -cm 1 Formatted $^)" == "x0" ] || exit 1

LINT_OPTS := -weak -fullinitblock
lint: src/*.c
	splint $(LINT_OPTS) $^

shellcheck.log: autogen
	shellcheck $^ | tee shellcheck.log
shellcheck: shellcheck.log
	[ "x$(shell wc -l $^)" == "x0 $^" ] || exit 1

commit:
	make shellcheck ctags check_style lint
	git ls-files -m
	[ "x$(shell git ls-files -m | wc -l)" == "x0" ] # no uncommitted files?
	make clean && CC=clang ./configure && make check distcheck
	make clean && ./configure && make check distcheck
	make prefix-installcheck
	make cov
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


EXTRA_DIST += $(TESTSUITE_AT) $(srcdir)/package.m4 $(TESTSUITE) atlocal atconfig
DISTCLEANFILES = atlocal atconfig
TESTSUITE_AT = testsuite.at
TESTSUITE = $(srcdir)/testsuite
test: check

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

LOCAL_INSTALL_PATH="$(PWD)/inst"
prefix-installcheck:
	-[ ! -f config.status.old ] && mv config.status config.status.old
	./configure --prefix=$(LOCAL_INSTALL_PATH)
	make install
	make installcheck
	@echo "---------------------"
	@echo "   installed files"
	@echo "---------------------"
	@find $(LOCAL_INSTALL_PATH) -type f
	@echo "---------------------"
	make uninstall
	@echo "all cleaned up?"
	find $(LOCAL_INSTALL_PATH) -type f
	[ "x$(shell find $(LOCAL_INSTALL_PATH) -type f | wc -l)" == "x0" ]
	mv config.status.old config.status

installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' MANPATH='$(mandir)' $(TESTSUITEFLAGS) -k 'install'
# don't run dependancy and unit tests

clean-local-check:
	test ! -f '$(TESTSUITE)' || $(SHELL) '$(TESTSUITE)' --clean
	-rm -f $(TESTSUITE)
	-rm -f $(srcdir)/package.m4

atconfig: $(top_builddir)/config.status
	$(top_builddir)/config.status $@

AUTOM4TE = $(SHELL) $(abs_srcdir)/build-aux/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(TESTSUITE_AT) $(srcdir)/package.m4
	$(AUTOTEST) -I '$(abs_srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@


# valgrind configuration
VALGRIND_BIN=valgrind
#VALGRIND_ENV=MPIWRAP_DEBUG=${MPIWRAP_DEBUG_ARGS} LD_PRELOAD=${MPI_PREFIX}/lib/valgrind/libmpiwrap-${MPI_PLATFORM}.so
VALGRIND_ARGS=--trace-children=yes --track-origins=yes --leak-check=full --show-leak-kinds=all
#--db-attach=yes
#MPI_BIN=mpirun
#MPI_ARGS=-n 4
VALGRIND=${VALGRIND_ENV} ${VALGRIND_BIN} ${VALGRIND_ARGS} ${MPI_BIN} ${MPI_ARGS}
valgrind:
	${VALGRIND} ./zedhat -t 9e-12 -f tests/testAXB-v73.mat

.PHONY: help
help:
	@echo "Useful make targets: make [target]"
	@echo "  <none> == all"
	@echo "  all              - build software (default)"
	@echo "  check            - run testsuite"
	@echo "  install          - install to --prefix=XXX"
	@echo "  installcheck     - test installed software"
	@echo "  uninstall        - uninstall from --prefix=XXX"
	@echo "  clean            - remove built files"
	@echo "For developers:"
	@echo "  distclean        - clean to distributed state"
	@echo "  maintainer-clean - delete all generated files, requires ./autogen to rebuild"
	@echo "  dist             - create a .tar.gz for distribution"
	@echo "  distcheck        - .tar, then untar in new dir and test"
	@echo "  prefix-installcheck - check for clean install/uninstall and run installcheck"
	@echo "  cov              - report test coverage"
	@echo "  ctags            - generate ctags for editing source files"
	@echo "  style            - run astyle on source files to tidy formatting"
	@echo "  lint             - run splint on source files"
