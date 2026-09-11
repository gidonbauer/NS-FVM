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

# = Stdpar =======================================
ifeq (${PARALLEL}, 1)
  ifeq (${COMP}, CLANG)
		CXX_FLAGS += -fexperimental-library -Wno-pass-failed
  else ifneq ($(filter GNU INTEL, ${COMP}),)  # Filter matches `GNU` and `INTEL`; true if not nothing is returned
		TBB_DIR ?= /opt/homebrew/opt/tbb
		TBB_INC = -I${TBB_DIR}/include
		TBB_LIB = -L${TBB_DIR}/lib -ltbb
		CXX_INC += ${TBB_INC}
		CXX_LIB += ${TBB_LIB}
  else ifeq (${COMP}, NVIDIA)
		CXX_FLAGS += -stdpar=gpu -Minfo=accel,par,stdpar -DIGOR_USE_CASSERT
  endif 
endif
# = Stdpar =======================================

# = HDF5 =========================================
HDF_DIR ?= /opt/homebrew/opt/hdf5
HDF_INC = -I${HDF_DIR}/include
HDF_LIB = -L${HDF_DIR}/lib -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
CXX_INC += ${HDF_INC}
CXX_LIB += ${HDF_LIB}
# = HDF5 =========================================
