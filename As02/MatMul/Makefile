# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -g -O1 -march=native
PAR_FLAG = -fopenmp
#GENERATE = -D GENERATE
GENERATE = 
TIMING = -D TIMING


APPS = matmul_basic matmul_opt matmul_omp matmul_cuda

.PHONY: all

all: $(APPS)

matmul_basic: matmul.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

matmul_opt: matmul_opt.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

matmul_omp: matmul_omp.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) -fopenmp $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

matmul_cuda: matmul_cuda.cu mmio.cpp MatrixParser/matrixparser.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^


.PHONY: clean
clean:
	rm -f $(APPS) *.o
