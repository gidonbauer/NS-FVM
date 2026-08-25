OBJ_DIR = obj
LIB_DIR = lib

OBJ_NAMES = ${OBJ_DIR}/precisions       \
            ${OBJ_DIR}/parameters       \
            ${OBJ_DIR}/customfftw3      \
            ${OBJ_DIR}/custompfft       \
            ${OBJ_DIR}/fft              \
            ${OBJ_DIR}/poisfft          \
            ${OBJ_DIR}/main             \
            ${OBJ_DIR}/c_binding        \
            ${OBJ_DIR}/f_mpi_comm_c2f   \
            ${OBJ_DIR}/c_binding_c

FFTW_DIR ?= /opt/homebrew/opt/fftw
# FFTW_INC = -I${FFTW_DIR}/include
FFTW_LIB = -L${FFTW_DIR}/lib -lfftw3 -lfftw3f
FFTW_LIB_OMP = -L${FFTW_DIR}/lib -lfftw3 -lfftw3f -lfftw3_omp # -lfftw3l_threads -lfftw3f_threads

all: ${LIB_DIR}/libpoisfft.dylib ${LIB_DIR}/libpoisfft_omp.dylib

${LIB_DIR}/libpoisfft.dylib: ${addsuffix .o, ${OBJ_NAMES}} | ${LIB_DIR}
	${FC} -dynamiclib -O3 -fbacktrace -fPIC $^ -Wl,-install_name,@rpath/libpoisfft.dylib -o $@ ${FFTW_LIB}

${LIB_DIR}/libpoisfft_omp.dylib: ${addsuffix _omp.o, ${OBJ_NAMES}} | ${LIB_DIR}
	${FC} -dynamiclib -O3 -fbacktrace -fopenmp -fPIC $^ -Wl,-install_name,@rpath/libpoisfft_omp.dylib -o $@ ${FFTW_LIB_OMP}

${OBJ_DIR}/%.o: src/%.f90 | ${OBJ_DIR}
	${FC} -c -cpp -g -O3 -fbacktrace -fPIC -I${OBJ_DIR} -J${OBJ_DIR} -o $@ $<

${OBJ_DIR}/%_omp.o: src/%.f90 | ${OBJ_DIR}
	${FC} -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I./${OBJ_DIR} -J./${OBJ_DIR} -o $@ $<

${OBJ_DIR}/%.o: src/%.c | ${OBJ_DIR}
	${CC} -O3 -c -fPIC -o $@ $<

${OBJ_DIR}/%_omp.o: src/%.c | ${OBJ_DIR}
	${CC} -O3 -c -fopenmp -fPIC -o $@ $<

${OBJ_DIR}/c_binding_c.o: src/c_binding.c | ${OBJ_DIR}
	${CC} -O3 -c -fPIC -o $@ $<

${OBJ_DIR}/c_binding_c_omp.o: src/c_binding.c | ${OBJ_DIR}
	${CC} -O3 -c -fopenmp -fPIC -o $@ $<

${OBJ_DIR}:
	mkdir -p $@

${LIB_DIR}:
	mkdir -p $@

clean:
	rm -fr ${OBJ_DIR} ${LIB_DIR}

.PHONY: all clean

# /opt/homebrew/bin/gfortran
# -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/lib/gcc/libpoisfft.dylib
# -dynamiclib
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/precisions.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/parameters.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/customfftw3.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/custompfft.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/fft.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/poisfft.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/main.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/f_mpi_comm_c2f.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding.os
# /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding_c.os

# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/precisions.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc precisions.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/parameters.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc parameters.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/customfftw3.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc customfftw3.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/custompfft.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc custompfft.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/fft.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc fft.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/poisfft.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc poisfft.f90
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/main.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc main.f90
# /opt/homebrew/opt/llvm/bin/clang -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/f_mpi_comm_c2f.os -c -fPIC f_mpi_comm_c2f.c
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding.os -c -cpp -g -O3 -fbacktrace -fopenmp -fPIC -I/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc -J/Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc c_binding.f90
# /opt/homebrew/opt/llvm/bin/clang -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding_c.os -c -fPIC c_binding.c
# /opt/homebrew/bin/gfortran -o /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/lib/gcc/libpoisfft.dylib -dynamiclib /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/precisions.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/parameters.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/customfftw3.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/custompfft.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/fft.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/poisfft.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/main.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/f_mpi_comm_c2f.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding.os /Users/gidonbauer/Documents/NS-FVM/Thirdparty/PoisFFT/bin/gcc/c_binding_c.os
