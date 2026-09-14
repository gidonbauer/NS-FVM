BENCH = ${addprefix bin/, ${basename ${wildcard bench/*.cpp}}}

bench: ${BENCH}

${BENCH}: bin/bench/%: bench/%.cpp ${HEADERS} | bin/bench
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin/bench:
	mkdir -p $@
