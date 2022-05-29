# See LICENSE.txt for license details.
CXX_FLAGS += -std=c++11 -O3 -Wall
PAR_FLAG = -fopenmp
#GENERATE = -D GENERATE
GENERATE = 
TIMING = -D TIMING

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif

APPS = matmul_cuda csr_cuda csc_cuda

.PHONY: all

all: $(APPS)

matmul_cuda: matmul_cuda.cu mmio.cpp MatrixParser/matrixparser.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

csr_cuda: csr.cu mmio.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

csc_cuda: csc.cu mmio.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

.PHONY: clean
clean:
	rm -f $(APPS)
