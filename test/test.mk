DEFAULT_BUILD_TESTS = test/Advection-Cartesian.cpp test/Advection-Polar.cpp test/Iterator.cpp test/Multigrid.cpp test/Polar-Couette.cpp test/Polar-Channel.cpp
DEFAULT_BUILD_TESTS := ${addprefix bin/, ${basename ${DEFAULT_BUILD_TESTS}}}

${DEFAULT_BUILD_TESTS}: bin/test/%: test/%.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-FFT: test/Taylor-Green.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-MG: test/Taylor-Green.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-FFT: test/Channel.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-MG: test/Channel.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

test/output bin/test:
	mkdir -p $@
