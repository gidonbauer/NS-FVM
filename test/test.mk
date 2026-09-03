TESTS = Taylor-Green-MG-8 Taylor-Green-MG-16 Taylor-Green-MG-64 \
        Taylor-Green-FFT-8 Taylor-Green-FFT-16 Taylor-Green-FFT-64 \
        Multigrid-32 Multigrid-64 Multigrid-128 Multigrid-512 Multigrid-1024 \
        Polar-Couette-8 Polar-Couette-16 Polar-Couette-32 \
        Channel-MG-16 Channel-MG-32 Channel-MG-64 \
        Channel-FFT-16 Channel-FFT-32 Channel-FFT-64

# ------------------------------------------------------------------------------
strip_digits = ${subst 9,,${subst 8,,${subst 7,,${subst 6,,${subst 5,,${subst 4,,${subst 3,,${subst 2,,${subst 1,,${subst 0,,$1}}}}}}}}}}
test_arg     = ${strip ${foreach l,${lastword ${subst -, ,$1}},${if ${call strip_digits,${l}},,${l}}}}
test_name    = ${patsubst %-${call test_arg,$1},%,$1}

TEST_NAMES = ${sort ${foreach t,${TESTS},${call test_name,${t}}}}
TEST_BINS  = ${addprefix bin/test/, ${TEST_NAMES}}

CUSTOM_TEST_BINS  = bin/test/Taylor-Green-FFT bin/test/Taylor-Green-MG bin/test/Channel-FFT bin/test/Channel-MG
GENERIC_TEST_BINS = ${filter-out ${CUSTOM_TEST_BINS}, ${TEST_BINS}}
# ------------------------------------------------------------------------------

test: ${addprefix test-, ${TESTS}}

.SECONDEXPANSION:
test-%: bin/test/$${call test_name,$$*}
	@printf "\033[32m[TEST]\033[0m Running test case $*...\n"
	@OMP_NUM_THREADS=4 $< ${call test_arg,$*} && printf "\033[32m[PASS]\033[0m $* finished successfully.\n\n" || printf "\033[31m[FAIL]\033[0m $* failed.\n\n"

${GENERIC_TEST_BINS}: bin/test/%: test/%.cpp ${HEADERS} | bin/test
	${CXX} ${CXX_FLAGS} ${CXX_INC} -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-FFT: test/Taylor-Green.cpp ${HEADERS} | bin/test
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Taylor-Green-MG: test/Taylor-Green.cpp ${HEADERS} | bin/test
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-FFT: test/Channel.cpp ${HEADERS} | bin/test
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DFFT_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test/Channel-MG: test/Channel.cpp ${HEADERS} | bin/test
	${CXX} ${CXX_FLAGS} ${CXX_INC} -DMG_POISSON=1 -o $@ $< ${CXX_LIB}

bin/test:
	mkdir -p $@

.PHONY: test
