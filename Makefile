TARGETS = ${addprefix bin/, ${basename ${wildcard *.cpp}}}
HEADERS = ${wildcard src/*.hpp}

CXX_FLAGS = -Wall -Wextra -pedantic -Wshadow -Wconversion -std=c++23
CXX_DEBUG_FLAGS = -O0 -g -fsanitize=address,undefined
CXX_RELEASE_FLAGS = -march=native -O3
CXX_FAST_FLAGS = ${CXX_RELEASE_FLAGS} -ffast-math -DNDEBUG -DIGOR_NDEBUG

DEBUG ?= 0
FAST ?= 0
PARALLEL ?= 0
ifeq (${DEBUG}, 1)
  CXX_FLAGS += ${CXX_DEBUG_FLAGS}
else ifeq (${FAST}, 1)
  CXX_FLAGS += ${CXX_FAST_FLAGS}
else
  CXX_FLAGS += ${CXX_RELEASE_FLAGS}
endif

ifeq (${PARALLEL}, 1)
	CXX_FLAGS += -DNS_FVM_PARALLEL
endif

CXX_INC := -I./src
CXX_LIB :=

IGOR_DIR ?= ${HOME}/opt/Igor
IGOR_INC = -I${IGOR_DIR}
CXX_INC += ${IGOR_INC}

POISFFT_DIR ?= ${HOME}/opt/PoisFFT
POISFFT_INC = -I${POISFFT_DIR}/src
POISFFT_LIB = -L${POISFFT_DIR}/lib/gcc/ -Wl,-rpath,${POISFFT_DIR}/lib/gcc/ -lpoisfft
CXX_INC += ${POISFFT_INC}
CXX_LIB += ${POISFFT_LIB}

ifeq (${PARALLEL}, 1)
	OMP_DIR ?= /opt/homebrew/opt/libomp/
	OMP_LIB = -L${OMP_DIR}/lib -lomp
	CXX_LIB += ${OMP_LIB}

	KOKKOS_DIR ?= /opt/homebrew/opt/kokkos
	KOKKOS_INC = -isystem${KOKKOS_DIR}/include
	KOKKOS_LIB = -L${KOKKOS_DIR}/lib -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkossimd
	CXX_FLAGS += -Xpreprocessor -fopenmp
	CXX_INC += ${KOKKOS_INC}
	CXX_LIB += ${KOKKOS_LIB}
endif

all: ${TARGETS}

bin/%: %.cpp ${HEADERS} | bin output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin:
	mkdir -p $@

output:
	mkdir -p $@

clean:
	rm -fr bin

.PHONY: all clean
