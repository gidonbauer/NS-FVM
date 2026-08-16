TARGETS = ${addprefix bin/, ${basename ${wildcard *.cpp}}}
HEADERS = ${wildcard src/*.hpp}

CXX_FLAGS = -Wall -Wextra -pedantic -Wshadow -Wconversion -std=c++23
CXX_DEBUG_FLAGS = -O0 -g -fsanitize=address,undefined
CXX_RELEASE_FLAGS = -march=native -O3
CXX_FAST_FLAGS = ${CXX_RELEASE_FLAGS} -ffast-math -DNDEBUG -DIGOR_NDEBUG

DEBUG ?= 0
FAST ?= 0
ifeq (${DEBUG}, 1)
  CXX_FLAGS += ${CXX_DEBUG_FLAGS}
else ifeq (${FAST}, 1)
  CXX_FLAGS += ${CXX_FAST_FLAGS}
else
  CXX_FLAGS += ${CXX_RELEASE_FLAGS}
endif

CXX_INC := -I./src

IGOR_DIR ?= ${HOME}/opt/Igor
IGOR_INC = -I${IGOR_DIR}
CXX_INC += ${IGOR_INC}

all: ${TARGETS}

bin/%: %.cpp ${HEADERS} | bin output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $<

bin:
	mkdir -p $@

output:
	mkdir -p $@

clean:
	rm -fr bin

.PHONY: all clean
