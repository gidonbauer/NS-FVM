BASENAME_CXX = ${notdir ${CXX}}

ifeq (${BASENAME_CXX}, clang++)
	COMP = CLANG
else ifneq ($(filter g++ g++-%, ${BASENAME_CXX}),)  # Filter matches all `g++` and `g++-*`; true if not nothing is returned
	COMP = GNU
else ifeq (${BASENAME_CXX}, icpx)
	COMP = INTEL
else ifeq (${BASENAME_CXX}, nvc++)
	COMP = NVIDIA
else
  ${error "Unknown C++ compiler `${CXX}`"}
endif

ifeq (${COMP}, CLANG)

	CXX_FLAGS           = -Wall -Wextra -pedantic -Wshadow -Wconversion -Winline -std=c++23
	CXX_RELEASE_FLAGS   = -march=native -mtune=native -O3
	CXX_FAST_FLAGS      = ${CXX_RELEASE_FLAGS} -ffast-math -DNDEBUG -DIGOR_NDEBUG
	CXX_DEBUG_FLAGS     = -O0 -g
	CXX_SANITIZER_FLAGS = -fsanitize=address,undefined

else ifeq (${COMP}, GNU)

	CXX_FLAGS           = -Wall -Wextra -pedantic -Wconversion -Wno-changes-meaning -std=c++23
	CXX_RELEASE_FLAGS   = -march=native -O3
	CXX_FAST_FLAGS      = ${CXX_RELEASE_FLAGS} -ffast-math -DNDEBUG -DIGOR_NDEBUG
	CXX_DEBUG_FLAGS     = -O0 -g
	CXX_SANITIZER_FLAGS = -fsanitize=address,undefined

else ifeq (${COMP}, INTEL)

	CXX_FLAGS           = -Wall -Wextra -pedantic -Wshadow -Wconversion -std=c++23
	CXX_RELEASE_FLAGS   = -O3 -march=native -mtune=native -fp-model precise
	CXX_FAST_FLAGS      = -O3 -march=native -mtune=native -fp-model fast=2 -ffast-math -DNDEBUG -DIGOR_NDEBUG
	CXX_DEBUG_FLAGS     = -O0 -g
	CXX_SANITIZER_FLAGS = -fsanitize=address,leak,undefined

else ifeq (${COMP}, NVIDIA)

	CXX_FLAGS           = -Wall -Wextra -pedantic -Wshadow -std=c++23
	CXX_RELEASE_FLAGS   = -O3 -fastsse -Mvect=simd:256,noassoc
	CXX_FAST_FLAGS      = -O3 -fast -fastsse -Mvect=simd:256 -DNDEBUG -DIGOR_NDEBUG
	CXX_DEBUG_FLAGS     = -O0 -g
	CXX_SANITIZER_FLAGS = -fsanitize=address,leak,undefined

else

  ${error "Unknown compiler `${COMP}`"}

endif

DEBUG    ?= 0
FAST     ?= 0
SANITIZE ?= 0
SCOREP   ?= 0
PARALLEL ?= 0

ifeq (${DEBUG}, 1)
  CXX_FLAGS += ${CXX_DEBUG_FLAGS}
else ifeq (${FAST}, 1)
  CXX_FLAGS += ${CXX_FAST_FLAGS}
else
  CXX_FLAGS += ${CXX_RELEASE_FLAGS}
endif

ifeq (${SANITIZE}, 1)
  CXX_FLAGS += ${CXX_SANITIZER_FLAGS}
endif

ifeq (${SCOREP}, 1)
  CXX_FLAGS += -g -fno-omit-frame-pointer
  CXX := scorep ${CXX}
endif

ifeq (${PARALLEL}, 1)
  CXX_FLAGS += -DNS_FVM_PARALLEL
endif
