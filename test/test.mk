bin/test/%: test/%.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-FFT: test/Taylor-Green.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-MG: test/Taylor-Green.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-FFT: test/Channel.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-MG: test/Channel.cpp ${HEADERS} | bin/test test/output
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

test/output:
	mkdir -p $@

bin/test:
	mkdir -p $@

.PHONY: test
