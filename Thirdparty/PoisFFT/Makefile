# ============================================================================
#  PoisFFT -- Makefile build system
#
#  Replaces the SCons build (src/SConstruct).  Builds four libraries:
#
#      lib/libpoisfft.a          lib/libpoisfft.(so|dylib)        serial
#      lib/libpoisfft_omp.a      lib/libpoisfft_omp.(so|dylib)    OpenMP
#
#  Supported Fortran compilers: GNU (gfortran), Intel (ifx/ifort),
#  NVIDIA (nvfortran) and LLVM Flang (flang / flang-new).
#
#  Quick start
#  -----------
#      make                      # auto-detect compiler, build everything
#      make COMP=intel           # force a compiler family
#      make FC=ftn COMP=cray     # explicit compiler, explicit family
#      make BUILD=fast           # relaxed-FP, aggressive optimisation
#      make BUILD=debug          # unoptimised, runtime checks
#      make NATIVE_ARCH=         # portable build, no CPU tuning
#      make test                 # build and run the bundled test programs
#      make install PREFIX=/opt/poisfft
#      make help                 # show all knobs
#
#  MPI/PFFT builds are not covered here: they additionally require the PFFT
#  library and its `pfft.f03' interface file, which is not shipped in src/.
#
#  Optimised builds tune for the build machine  (IMPORTANT)
#  --------------------------------------------------------
#  Both optimised build types -- the default BUILD=release and BUILD=fast --
#  compile for the CPU they are built on: -mcpu=native on AArch64,
#  -march=native elsewhere, -xHost for Intel, -tp=native for NVIDIA.
#
#  On a cluster whose login node differs from its compute nodes that produces
#  binaries which die with an illegal instruction.  Either build on a compute
#  node, or name the target explicitly:
#
#      make NATIVE_ARCH=-march=znver3      # build for a specific CPU
#      make NATIVE_ARCH=                   # portable, no CPU tuning
#
#  NATIVE_ARCH is the single knob for this and applies to release and fast.
#  BUILD=debug never adds it.
#
#  A note on BUILD=fast
#  --------------------
#  On top of the CPU tuning it relaxes IEEE floating point (-ffast-math /
#  -fp-model=fast / -Mfprelaxed / -hfp3).  Besides costing a little accuracy,
#  GCC links crtfastmath.o into the resulting library, which switches
#  flush-to-zero on for the whole process -- including the calling program.
#  Stay on BUILD=release if the caller depends on denormals or strict IEEE
#  semantics.
#
#  Its flags are spelled out rather than using -Ofast, on purpose: clang/flang
#  deprecated -Ofast in LLVM 22, and gfortran's -Ofast quietly implies
#  -fstack-arrays, which moves array temporaries onto the (small, per-thread)
#  OpenMP stacks and can overflow them on large grids.
#
#  Link-time optimisation is deliberately *not* part of it: -flto needs a
#  plugin-aware archiver for the static libraries (AR=gcc-ar / xiar) and is
#  fragile across the Fortran/C/C++ boundary.  Add it deliberately if you want
#  it: `make BUILD=fast EXTRA_FFLAGS=-flto EXTRA_CFLAGS=-flto AR=gcc-ar'.
# ============================================================================

MAKEFLAGS += --no-builtin-rules --no-builtin-variables
.SUFFIXES:

SRC_DIR   := src
BUILD_DIR := build
LIB_DIR   := lib
BIN_DIR   := bin

BUILD     ?= release
PREFIX    ?= /usr/local

# ----------------------------------------------------------------------------
# Host platform
# ----------------------------------------------------------------------------
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Flag that tunes for the machine we are building on, used by the optimised
# build types (release and fast).
# GCC/LLVM spell it -mcpu=native on AArch64 and -march=native elsewhere;
# Intel and NVIDIA have their own spelling (set in their family blocks).
ifneq (,$(filter arm64 aarch64,$(UNAME_M)))
  HOST_ARCH_FLAG := -mcpu=native
else
  HOST_ARCH_FLAG := -march=native
endif

# `-shared' is understood by every supported driver, including on Darwin,
# where the GCC/LLVM drivers translate it to -dynamiclib.  flang rejects
# -dynamiclib outright, so do not use it.
ifeq ($(UNAME_S),Darwin)
  SHLIB_EXT   := dylib
  SHLIB_FLAG   = -shared -Wl,-install_name,@rpath/$(notdir $@)
else
  SHLIB_EXT   := so
  SHLIB_FLAG   = -shared -Wl,-soname,$(notdir $@)
endif

# ----------------------------------------------------------------------------
# Compiler selection
#
# `make' pre-defines FC/CC/CXX to f77/cc/g++; ignore those defaults so that we
# can auto-detect, but honour anything the user or the environment sets.
# ----------------------------------------------------------------------------
ifeq ($(origin FC),default)
  FC :=
endif
ifeq ($(origin CC),default)
  CC :=
endif
ifeq ($(origin CXX),default)
  CXX :=
endif

which = $(shell command -v $(1) 2>/dev/null)
first_of = $(firstword $(foreach c,$(1),$(if $(call which,$(c)),$(c))))

# ----------------------------------------------------------------------------
# Compiler family detection (override with COMP=gnu|intel|nvidia|flang|cray)
#
# When COMP is given explicitly the family block below picks that family's
# default compiler, so only probe the PATH when we have to guess.
# ----------------------------------------------------------------------------
COMP ?= auto

ifeq ($(COMP),auto)
  ifeq ($(strip $(FC)),)
    FC := $(call first_of,gfortran ifx nvfortran flang flang-new ifort ftn)
  endif
  ifeq ($(strip $(FC)),)
    $(error No Fortran compiler found in PATH. Set FC=<compiler> and COMP=<family>)
  endif

  FC_BANNER := $(shell $(FC) --version 2>&1 | head -n 2 | tr '\n' ' ')
  ifneq (,$(findstring GNU Fortran,$(FC_BANNER)))
    COMP := gnu
  else ifneq (,$(findstring flang,$(FC_BANNER)))
    COMP := flang
  else ifneq (,$(findstring IFX,$(FC_BANNER)))
    COMP := intel
  else ifneq (,$(findstring IFORT,$(FC_BANNER)))
    COMP := intel
  else ifneq (,$(findstring ifort,$(FC_BANNER)))
    COMP := intel
  else ifneq (,$(findstring nvfortran,$(FC_BANNER)))
    COMP := nvidia
  else ifneq (,$(findstring NVIDIA,$(FC_BANNER)))
    COMP := nvidia
  else ifneq (,$(findstring pgf,$(FC_BANNER)))
    COMP := nvidia
  else ifneq (,$(findstring Cray,$(FC_BANNER)))
    COMP := cray
  else
    $(warning Could not identify the Fortran compiler '$(FC)' from its version banner)
    $(warning Assuming GNU-style flags; override with COMP=gnu|intel|nvidia|flang|cray)
    COMP := gnu
  endif
endif

# ----------------------------------------------------------------------------
# Per-family flags
#
#   FC_PP      preprocess Fortran sources with cpp   (sources use #include/#if)
#   FC_PIC     position independent code
#   FC_OMP     enable OpenMP
#   FC_MODOUT  where to write .mod files  (uses $(@D), i.e. the object dir)
#   NATIVE_ARCH tune for this machine's CPU; part of every optimised build
#   FC_OPT/CC_OPT   flags for BUILD=release
#   FC_FAST/CC_FAST flags for BUILD=fast
#   FC_DBG     flags for BUILD=debug
#   FC_WARN    warnings, only added with WARNINGS=1
# ----------------------------------------------------------------------------
ifeq ($(COMP),gnu)
  FC        := $(if $(FC),$(FC),gfortran)
  CC        := $(if $(CC),$(CC),$(call first_of,gcc clang cc))
  CXX       := $(if $(CXX),$(CXX),$(call first_of,g++ clang++ c++))
  FC_PP     := -cpp
  FC_PIC    := -fPIC
  FC_OMP    := -fopenmp
  FC_MODOUT  = -J$(@D)
  NATIVE_ARCH ?= $(HOST_ARCH_FLAG)
  FC_OPT    := -O3 $(NATIVE_ARCH) -g -fbacktrace
  CC_OPT    := -O3 $(NATIVE_ARCH) -g
  FC_FAST   := -O3 -ffast-math -fno-protect-parens $(NATIVE_ARCH) -funroll-loops -g -fbacktrace
  CC_FAST   := -O3 -ffast-math $(NATIVE_ARCH) -funroll-loops -g
  FC_DBG    := -O0 -g -fbacktrace -fcheck=bounds,pointer -ffpe-trap=zero,overflow,invalid
  FC_WARN   := -Wall -Wsurprising -Wno-unused-function -Wno-unused-dummy-argument
  FC_SANIT  := -fsanitize=undefined,address

else ifeq ($(COMP),intel)
  FC        := $(if $(FC),$(FC),ifx)
  CC        := $(if $(CC),$(CC),$(call first_of,icx icc cc))
  CXX       := $(if $(CXX),$(CXX),$(call first_of,icpx icpc c++))
  FC_PP     := -fpp
  FC_PIC    := -fPIC
  FC_OMP    := -qopenmp
  FC_MODOUT  = -module $(@D)
  NATIVE_ARCH ?= -xHost
  FC_OPT    := -O3 $(NATIVE_ARCH) -g -traceback
  CC_OPT    := -O3 $(NATIVE_ARCH) -g
  FC_FAST   := -O3 -fp-model=fast $(NATIVE_ARCH) -funroll-loops -g -traceback
  CC_FAST   := -O3 -fp-model=fast $(NATIVE_ARCH) -funroll-loops -g
  FC_DBG    := -O0 -g -traceback -check all,noarg_temp_created
  FC_WARN   := -warn all -warn nounused
  FC_SANIT  := -fsanitize=address

else ifeq ($(COMP),nvidia)
  FC        := $(if $(FC),$(FC),nvfortran)
  CC        := $(if $(CC),$(CC),$(call first_of,nvc cc gcc))
  CXX       := $(if $(CXX),$(CXX),$(call first_of,nvc++ c++ g++))
  FC_PP     := -Mpreprocess
  FC_PIC    := -fPIC
  FC_OMP    := -mp
  FC_MODOUT  = -module $(@D)
  NATIVE_ARCH ?= -tp=native
  FC_OPT    := -O3 $(NATIVE_ARCH) -g -traceback
  CC_OPT    := -O3 $(NATIVE_ARCH) -g
  FC_FAST   := -fast -O3 $(NATIVE_ARCH) -Mvect=simd -Mfprelaxed -g -traceback
  CC_FAST   := -fast -O3 $(NATIVE_ARCH) -Mvect=simd -Mfprelaxed -g
  FC_DBG    := -O0 -g -traceback -Mbounds -Mchkptr
  FC_WARN   := -Minform=warn
  FC_SANIT  :=

else ifeq ($(COMP),flang)
  FC        := $(if $(FC),$(FC),flang)
  CC        := $(if $(CC),$(CC),$(call first_of,clang cc gcc))
  CXX       := $(if $(CXX),$(CXX),$(call first_of,clang++ c++ g++))
  FC_PP     := -cpp
  FC_PIC    := -fPIC
  FC_OMP    := -fopenmp
  FC_MODOUT  = -J$(@D)
  NATIVE_ARCH ?= $(HOST_ARCH_FLAG)
  FC_OPT    := -O3 $(NATIVE_ARCH) -g
  CC_OPT    := -O3 $(NATIVE_ARCH) -g
  FC_FAST   := -O3 -ffast-math $(NATIVE_ARCH) -funroll-loops -g
  CC_FAST   := -O3 -ffast-math $(NATIVE_ARCH) -funroll-loops -g
  FC_DBG    := -O0 -g
  FC_WARN   := -Wall
  FC_SANIT  := -fsanitize=address

else ifeq ($(COMP),cray)
  # Cray compiler wrappers (ftn/cc/CC); the wrapper picks the real compiler.
  FC        := $(if $(FC),$(FC),ftn)
  CC        := $(if $(CC),$(CC),cc)
  CXX       := $(if $(CXX),$(CXX),CC)
  FC_PP     := -eZ
  FC_PIC    := -fPIC
  FC_OMP    := -fopenmp
  FC_MODOUT  = -J$(@D)
  NATIVE_ARCH ?=
  FC_OPT    := -O2 $(NATIVE_ARCH) -g
  CC_OPT    := -O2 $(NATIVE_ARCH) -g
  FC_FAST   := -O3 -hfp3 $(NATIVE_ARCH) -g
  CC_FAST   := -O3 -ffast-math $(NATIVE_ARCH) -g
  FC_DBG    := -O0 -g -Rb
  FC_WARN   :=
  FC_SANIT  :=

else
  $(error Unknown COMP='$(COMP)'. Use one of: gnu intel nvidia flang cray)
endif

CC  := $(if $(strip $(CC)),$(CC),cc)
CXX := $(if $(strip $(CXX)),$(CXX),c++)

ifeq ($(strip $(FC)),)
  $(error No Fortran compiler found. Set FC=<compiler> explicitly)
endif

# Extra flags needed to *link* an OpenMP build.  Normally the compiler driver
# handles this itself; Homebrew's flang is a known exception -- it accepts
# -fopenmp but does not tell the linker where libomp lives.
ifeq ($(COMP)-$(UNAME_S),flang-Darwin)
  LIBOMP_DIR ?= $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(strip $(LIBOMP_DIR)),)
    OMP_LDFLAGS ?= -L$(LIBOMP_DIR)/lib -Wl,-rpath,$(LIBOMP_DIR)/lib
  endif
endif
OMP_LDFLAGS ?=

# ----------------------------------------------------------------------------
# Build type
# ----------------------------------------------------------------------------
# FAST=1 was an earlier spelling of BUILD=fast.  Refuse it rather than
# silently handing back a plain release build.
ifdef FAST
  $(error FAST is no longer a separate switch; use BUILD=fast)
endif

ifeq ($(BUILD),release)
  FC_BUILDFLAGS := $(FC_OPT)
  CC_BUILDFLAGS := $(CC_OPT)
else ifeq ($(BUILD),fast)
  FC_BUILDFLAGS := $(FC_FAST)
  CC_BUILDFLAGS := $(CC_FAST)
else ifeq ($(BUILD),debug)
  FC_BUILDFLAGS := $(FC_DBG)
  CC_BUILDFLAGS := -O0 -g
else
  $(error Unknown BUILD='$(BUILD)'. Use BUILD=release, BUILD=fast or BUILD=debug)
endif

ifeq ($(WARNINGS),1)
  FC_BUILDFLAGS += $(FC_WARN)
  CC_BUILDFLAGS += -Wall -Wextra
endif

ifeq ($(SANITIZE),1)
  ifeq ($(strip $(FC_SANIT)),)
    $(warning SANITIZE=1 is not supported for COMP=$(COMP); ignored)
  else
    FC_BUILDFLAGS += $(FC_SANIT)
    CC_BUILDFLAGS += $(FC_SANIT)
    LD_SANIT      := $(FC_SANIT)
  endif
endif

# ----------------------------------------------------------------------------
# FFTW
#
# PoisFFT declares the FFTW interfaces itself (src/customfftw3.f90), so no FFTW
# headers are needed -- only the libraries, and only when linking.
# ----------------------------------------------------------------------------
comma := ,

FFTW_DIR ?=
ifneq ($(strip $(FFTW_DIR)),)
  FFTW_LDFLAGS := -L$(FFTW_DIR)/lib
else
  FFTW_LDFLAGS := $(shell pkg-config --libs-only-L fftw3 2>/dev/null)
endif
# Bake the FFTW location into the shared library so it is found at run time.
FFTW_RPATH ?= $(patsubst -L%,-Wl$(comma)-rpath$(comma)%,$(FFTW_LDFLAGS))

FFTW_LIBS     ?= -lfftw3 -lfftw3f
FFTW_OMP_LIBS ?= -lfftw3_omp

LINK_SERIAL := $(FFTW_LDFLAGS) $(FFTW_RPATH) $(FFTW_LIBS)
LINK_OMP    := $(FFTW_LDFLAGS) $(FFTW_RPATH) $(FFTW_LIBS) $(FFTW_OMP_LIBS) $(OMP_LDFLAGS)

# ----------------------------------------------------------------------------
# Assembled flags
# ----------------------------------------------------------------------------
FFLAGS  = $(FC_PP) $(FC_PIC) $(FC_BUILDFLAGS) $(FC_MODOUT) -I$(@D) -I$(SRC_DIR) $(EXTRA_FFLAGS)
CFLAGS  = -fPIC $(CC_BUILDFLAGS) -I$(SRC_DIR) $(EXTRA_CFLAGS)
LDFLAGS = $(LD_SANIT) $(EXTRA_LDFLAGS)

AR         ?= ar
POISFFT_ARFLAGS ?= rcs

# ----------------------------------------------------------------------------
# Objects
#
# Object order matters for the archive/link: modules must be compiled in
# dependency order (enforced by the explicit prerequisites further down).
#
# f_mpi_comm_c2f.c and missing_pfft.c are MPI/PFFT-only and compile to nothing
# without -DMPI, so they are left out.
# ----------------------------------------------------------------------------
F_OBJ_NAMES := precisions parameters customfftw3 custompfft fft poisfft main c_binding
C_OBJ_NAMES := c_binding_c
OBJ_NAMES   := $(F_OBJ_NAMES) $(C_OBJ_NAMES)

SER_DIR := $(BUILD_DIR)/$(COMP)/serial
OMP_DIR := $(BUILD_DIR)/$(COMP)/omp

SER_OBJS := $(addprefix $(SER_DIR)/,$(addsuffix .o,$(OBJ_NAMES)))
OMP_OBJS := $(addprefix $(OMP_DIR)/,$(addsuffix .o,$(OBJ_NAMES)))

LIB_STATIC     := $(LIB_DIR)/libpoisfft.a
LIB_STATIC_OMP := $(LIB_DIR)/libpoisfft_omp.a
LIB_SHARED     := $(LIB_DIR)/libpoisfft.$(SHLIB_EXT)
LIB_SHARED_OMP := $(LIB_DIR)/libpoisfft_omp.$(SHLIB_EXT)

STATIC_LIBS := $(LIB_STATIC) $(LIB_STATIC_OMP)
SHARED_LIBS := $(LIB_SHARED) $(LIB_SHARED_OMP)

# ============================================================================
# Targets
# ============================================================================
.PHONY: all static shared serial omp test check clean distclean install help info
.DEFAULT_GOAL := all

all: static shared

static: $(STATIC_LIBS)
shared: $(SHARED_LIBS)

serial: $(LIB_STATIC) $(LIB_SHARED)
omp:    $(LIB_STATIC_OMP) $(LIB_SHARED_OMP)

# --- libraries --------------------------------------------------------------
$(LIB_STATIC): $(SER_OBJS) | $(LIB_DIR)
	$(AR) $(POISFFT_ARFLAGS) $@ $^

$(LIB_STATIC_OMP): $(OMP_OBJS) | $(LIB_DIR)
	$(AR) $(POISFFT_ARFLAGS) $@ $^

$(LIB_SHARED): $(SER_OBJS) | $(LIB_DIR)
	$(FC) $(SHLIB_FLAG) $(FC_PIC) $(LDFLAGS) $^ -o $@ $(LINK_SERIAL)

$(LIB_SHARED_OMP): $(OMP_OBJS) | $(LIB_DIR)
	$(FC) $(SHLIB_FLAG) $(FC_PIC) $(FC_OMP) $(LDFLAGS) $^ -o $@ $(LINK_OMP)

# --- compilation ------------------------------------------------------------
$(SER_DIR)/%.o: $(SRC_DIR)/%.f90 | $(SER_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OMP_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OMP_DIR)
	$(FC) $(FFLAGS) $(FC_OMP) -c $< -o $@

$(SER_DIR)/%.o: $(SRC_DIR)/%.c | $(SER_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OMP_DIR)/%.o: $(SRC_DIR)/%.c | $(OMP_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# c_binding.c and c_binding.f90 would collide on the object name c_binding.o.
$(SER_DIR)/c_binding_c.o: $(SRC_DIR)/c_binding.c | $(SER_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OMP_DIR)/c_binding_c.o: $(SRC_DIR)/c_binding.c | $(OMP_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# --- inter-module dependencies ---------------------------------------------
# Fortran modules must be built in order, and the -inc.f90 files are pulled in
# textually by the preprocessor, so they are prerequisites too.
define POISFFT_DEPS
$(1)/custompfft.o: $(1)/customfftw3.o
$(1)/fft.o:       $(1)/precisions.o $(1)/customfftw3.o $(1)/custompfft.o \
                  $(SRC_DIR)/fft-inc.f90 $(SRC_DIR)/plan_new-inc.f90 \
                  $(SRC_DIR)/allocate_fftw-inc.f90
$(1)/poisfft.o:   $(1)/precisions.o $(1)/parameters.o $(1)/fft.o \
                  $(SRC_DIR)/poisfft-inc.f90 $(SRC_DIR)/poisfft-solvers-inc.f90 \
                  $(SRC_DIR)/poisfft-solvers-nonuniform_z-inc.f90 \
                  $(SRC_DIR)/loop_nest_3d.f90
$(1)/main.o:      $(1)/precisions.o $(1)/parameters.o $(1)/poisfft.o $(1)/custompfft.o
$(1)/c_binding.o: $(1)/main.o \
                  $(SRC_DIR)/c_new-inc.f90 $(SRC_DIR)/c_execute-inc.f90 \
                  $(SRC_DIR)/c_finalize-inc.f90
$(1)/c_binding_c.o: $(SRC_DIR)/poisfft.h
endef

$(eval $(call POISFFT_DEPS,$(SER_DIR)))
$(eval $(call POISFFT_DEPS,$(OMP_DIR)))

# --- directories ------------------------------------------------------------
$(SER_DIR) $(OMP_DIR) $(LIB_DIR) $(BIN_DIR):
	@mkdir -p $@

# ============================================================================
# Test programs (from the former SConstruct `test' alias)
#
# The Fortran test links the static library (the Fortran driver supplies its
# own runtime); the C/C++ tests link the shared library, which already records
# its dependency on the Fortran runtime.
# ============================================================================
TESTS := $(BIN_DIR)/testpoisson $(BIN_DIR)/testpoisson_omp \
         $(BIN_DIR)/c_testpoisson $(BIN_DIR)/cc_testpoisson

test check: $(TESTS)
	@for t in $(TESTS); do \
	  echo "=== $$t ==="; \
	  DYLD_LIBRARY_PATH=$(CURDIR)/$(LIB_DIR) LD_LIBRARY_PATH=$(CURDIR)/$(LIB_DIR) $$t || exit 1; \
	done
	@echo "=== all tests passed ==="

$(BIN_DIR)/testpoisson: $(SRC_DIR)/testpoisson.f90 $(LIB_STATIC) | $(BIN_DIR)
	$(FC) $(FC_PP) $(FC_BUILDFLAGS) -I$(SER_DIR) $(LDFLAGS) $< $(LIB_STATIC) -o $@ $(LINK_SERIAL)

$(BIN_DIR)/testpoisson_omp: $(SRC_DIR)/testpoisson.f90 $(LIB_STATIC_OMP) | $(BIN_DIR)
	$(FC) $(FC_PP) $(FC_BUILDFLAGS) $(FC_OMP) -I$(OMP_DIR) $(LDFLAGS) $< $(LIB_STATIC_OMP) -o $@ $(LINK_OMP)

$(BIN_DIR)/c_testpoisson: $(SRC_DIR)/testpoisson.c $(LIB_SHARED) | $(BIN_DIR)
	$(CC) $(CC_BUILDFLAGS) -I$(SRC_DIR) $(LDFLAGS) $< -o $@ \
	    -L$(LIB_DIR) -Wl,-rpath,$(CURDIR)/$(LIB_DIR) -lpoisfft -lm

$(BIN_DIR)/cc_testpoisson: $(SRC_DIR)/testpoisson.cc $(LIB_SHARED) | $(BIN_DIR)
	$(CXX) $(CC_BUILDFLAGS) -I$(SRC_DIR) $(LDFLAGS) $< -o $@ \
	    -L$(LIB_DIR) -Wl,-rpath,$(CURDIR)/$(LIB_DIR) -lpoisfft -lm

# ============================================================================
# install / clean / help
# ============================================================================
install: all
	@mkdir -p $(DESTDIR)$(PREFIX)/lib
	@mkdir -p $(DESTDIR)$(PREFIX)/include/poisfft/serial
	@mkdir -p $(DESTDIR)$(PREFIX)/include/poisfft/omp
	cp $(STATIC_LIBS) $(SHARED_LIBS) $(DESTDIR)$(PREFIX)/lib/
	cp $(SRC_DIR)/poisfft.h $(DESTDIR)$(PREFIX)/include/
	cp $(SER_DIR)/*.mod $(DESTDIR)$(PREFIX)/include/poisfft/serial/
	cp $(OMP_DIR)/*.mod $(DESTDIR)$(PREFIX)/include/poisfft/omp/
	@echo "Installed into $(DESTDIR)$(PREFIX)"

clean:
	rm -rf $(BUILD_DIR) $(LIB_DIR) $(BIN_DIR)

distclean: clean
	rm -rf obj $(SRC_DIR)/.sconsign.dblite

info:
	@echo "COMP        = $(COMP)"
	@echo "FC          = $(FC)   ($(shell $(FC) --version 2>&1 | head -n 1))"
	@echo "CC          = $(CC)"
	@echo "CXX         = $(CXX)"
	@echo "BUILD       = $(BUILD)"
	@echo "FFLAGS      = $(FC_PP) $(FC_PIC) $(FC_BUILDFLAGS) <moddir> -I$(SRC_DIR) $(EXTRA_FFLAGS)"
	@echo "OpenMP flag = $(FC_OMP)"
	@echo "FFTW link   = $(LINK_OMP)"
	@echo "objects     -> $(SER_DIR), $(OMP_DIR)"
	@echo "libraries   -> $(STATIC_LIBS) $(SHARED_LIBS)"

help:
	@echo "PoisFFT build system"
	@echo
	@echo "Targets:"
	@echo "  all        static + shared, serial + OpenMP (default)"
	@echo "  static     lib/libpoisfft.a and lib/libpoisfft_omp.a"
	@echo "  shared     lib/libpoisfft.$(SHLIB_EXT) and lib/libpoisfft_omp.$(SHLIB_EXT)"
	@echo "  serial     both flavours of the serial library only"
	@echo "  omp        both flavours of the OpenMP library only"
	@echo "  test       build and run the bundled test programs"
	@echo "  install    install to \$$PREFIX (default $(PREFIX))"
	@echo "  info       show the detected toolchain and flags"
	@echo "  clean      remove build/, lib/ and bin/"
	@echo "  distclean  clean, plus leftovers from the old SCons build"
	@echo
	@echo "Variables:"
	@echo "  COMP=gnu|intel|nvidia|flang|cray   compiler family (auto-detected: $(COMP))"
	@echo "  FC=, CC=, CXX=                     explicit compilers"
	@echo "  BUILD=release|fast|debug           default: release"
	@echo "    release   optimised, IEEE conformant"
	@echo "              ($(COMP): $(FC_OPT))"
	@echo "    fast      aggressive, relaxed floating point"
	@echo "              ($(COMP): $(FC_FAST))"
	@echo "    debug     unoptimised, runtime checks"
	@echo "              ($(COMP): $(FC_DBG))"
	@echo "  NATIVE_ARCH=<flag>                 CPU tuning for every optimised build,"
	@echo "                                     default '$(NATIVE_ARCH)'.  Set it explicitly"
	@echo "                                     (or empty) when the build machine is not"
	@echo "                                     the run machine, e.g. an HPC login node."
	@echo "  WARNINGS=1                         enable the warning flag set"
	@echo "  SANITIZE=1                         enable ASan/UBSan where supported"
	@echo "  FFTW_DIR=<prefix>                  FFTW install prefix (else pkg-config)"
	@echo "  FFTW_LIBS=, FFTW_OMP_LIBS=         override the FFTW link line"
	@echo "  OMP_LDFLAGS=                       extra flags to link an OpenMP build"
	@echo "  AR=, POISFFT_ARFLAGS=              archiver (default: $(AR) $(POISFFT_ARFLAGS))"
	@echo "  EXTRA_FFLAGS=, EXTRA_CFLAGS=, EXTRA_LDFLAGS="
	@echo "  PREFIX=, DESTDIR=                  install location"
	@echo
	@echo "Linking against PoisFFT:"
	@echo "  shared, serial:  -L$(LIB_DIR) -Wl,-rpath,\$$(pwd)/$(LIB_DIR) -lpoisfft"
	@echo "  shared, OpenMP:  -L$(LIB_DIR) -Wl,-rpath,\$$(pwd)/$(LIB_DIR) -lpoisfft_omp"
	@echo "  static, serial:  $(LIB_STATIC) $(LINK_SERIAL)"
	@echo "  static, OpenMP:  $(LIB_STATIC_OMP) $(LINK_OMP) $(FC_OMP)"
	@echo "  C/C++ callers additionally need -I$(SRC_DIR) for poisfft.h, and when"
	@echo "  linking a static build with a C/C++ driver, the Fortran runtime."
