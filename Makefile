TARGETS = ${addprefix bin/, ${basename ${wildcard *.cpp}}}
POISFFT_TARGETS = bin/advection-diffusion bin/mac bin/lid-driven-cavity bin/taylor-green
NON_POISFFT_TARGETS = ${filter-out ${POISFFT_TARGETS}, ${TARGETS}}
HEADERS = ${wildcard src/*.hpp}

include Makefiles/compiler_flags.mk
include Makefiles/libs.mk

all: ${TARGETS}

${NON_POISFFT_TARGETS}: bin/%: %.cpp ${HEADERS} | bin output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

${POISFFT_TARGETS}: bin/%: %.cpp ${HEADERS} | bin output
	${CXX} ${CXX_FLAGS} ${CXX_INC} ${POISFFT_INC} -o $@ $< ${CXX_LIB} ${POISFFT_LIB}

bin output: %:
	mkdir -p $@

clean:
	rm -fr bin

include test/test.mk
include bench/bench.mk

.PHONY: all clean
