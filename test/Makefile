## Directories:
# Package source files directory.
SRCDIR       := ../src
# Binary files directory.
BINDIR       := ./bin
# Make the binary files directory.
$(shell mkdir -p ${BINDIR} >/dev/null)
# Dependencies Makefiles' directory.
DPMDIR       := .d
# Make the dependencies Makefiles' directory.
$(shell mkdir -p ${DPMDIR} >/dev/null)

#Makefile absolute path.
TEST_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

##  Build variables:
# Compiler.
CXX          ?= g++
# Compiler flags.
MACROS       ?= -DTESTDIR=\"${TEST_PATH}\"
CPPFLAGS      = -Wall -pedantic -std=c++14 ${MACROS} ${LIBS_CFLAGS}
CXXFLAGS     ?= -g
LIBS_CFLAGS  := -I${SRCDIR}/ \
                -Iinclude/ \
                $(shell pkg-config --cflags seqan-2)
LDFLAG       := -lz \
                -lbz2 \
		-lpthread

##  Specifying source files to automatically generate dependencies:
# Consider all *.cc files in the src directory as 'source' file.
SOURCES      := $(wildcard *.cc)

##  Targets:
TARGETS       := ${BINDIR}/kseq++_test ${BINDIR}/kseq++_bench ${BINDIR}/seqio_test

##  Recipes:
COMPILE.cxx   = ${CXX} -c ${CXXFLAGS} ${CPPFLAGS} -o $@
LINK          = ${CXX} $^ ${LDFLAG} -o $@
MAKE_DMAKES   = ${CXX} ${CPPFLAGS} -MM -MG -MF $@

# Specifying phony targets.
.PHONY: all clean distclean

# building target.
all: ${TARGETS}

${BINDIR}/kseq++_test: kseq++_test.o
	${LINK}

${BINDIR}/kseq++_bench: kseq++_bench.o
	${LINK}

${BINDIR}/seqio_test: seqio_test.o
	${LINK}

# general rule for compiling cc files.
%.o: %.cc
	${COMPILE.cxx} $<

clean:
	rm -f *.o *~

distclean: clean
	rm -f -r ${BINDIR}; \
	rm -f -r ${DPMDIR}

# generating sources' prerequisites automatically by remaking *.d Makefiles.
${DPMDIR}/%.d: %.cc
	${MAKE_DMAKES} -MT $*.o -MT $@ $<

# include *.d Makefiles.
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.d})
endif
endif
