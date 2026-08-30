TARGETS = ${addprefix bin/, ${basename ${wildcard *.cpp}}}
HEADERS = ${wildcard src/*.hpp}

include Makefiles/compiler_flags.mk
include Makefiles/libs.mk

all: ${TARGETS}

${TARGETS}: bin/%: %.cpp ${HEADERS} | bin output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin:
	mkdir -p $@

output:
	mkdir -p $@

clean:
	rm -fr bin

include test/test.mk

.PHONY: all clean
