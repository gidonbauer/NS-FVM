CXX_INC := -I./src
CXX_LIB :=

# = Igor =========================================
IGOR_DIR ?= ${HOME}/opt/Igor
IGOR_INC = -I${IGOR_DIR}
CXX_INC += ${IGOR_INC}
# = Igor =========================================

# = cpptrace =====================================
ifeq (${BACKTRACE}, 1)
	CPPTRACE_DIR ?= /opt/homebrew/opt/cpptrace/
	CPPTRACE_INC = -I${CPPTRACE_DIR}/include -DIGOR_USE_CPPTRACE -g
	CPPTRACE_LIB = -L${CPPTRACE_DIR}/lib -Wl,-rpath,${CPPTRACE_DIR}/lib -lcpptrace
	CXX_INC += ${CPPTRACE_INC}
	CXX_LIB += ${CPPTRACE_LIB}
endif
# = cpptrace =====================================

# = PoisFFT ======================================
POISFFT_DIR ?= ./Thirdparty/PoisFFT
POISFFT_INC = -I${POISFFT_DIR}/src
POISFFT_LIB = -L${POISFFT_DIR}/lib/ -Wl,-rpath,${POISFFT_DIR}/lib/
ifeq (${PARALLEL}, 1)
	# POISFFT_LIB += -lpoisfft_omp
	POISFFT_LIB += -lpoisfft
else
	POISFFT_LIB += -lpoisfft
endif

CXX_INC += ${POISFFT_INC}
CXX_LIB += ${POISFFT_LIB}
# = PoisFFT ======================================

# = Kokkos =======================================
ifeq (${PARALLEL}, 1)
	OMP_DIR ?= /opt/homebrew/opt/libomp
	OMP_LIB = -L${OMP_DIR}/lib -lomp
	CXX_LIB += ${OMP_LIB}

	KOKKOS_DIR ?= /opt/homebrew/opt/kokkos
	KOKKOS_INC = -isystem${KOKKOS_DIR}/include
	KOKKOS_LIB = -L${KOKKOS_DIR}/lib -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkossimd
	CXX_FLAGS += -Xpreprocessor -fopenmp
	CXX_INC += ${KOKKOS_INC}
	CXX_LIB += ${KOKKOS_LIB}
endif
# = Kokkos =======================================
