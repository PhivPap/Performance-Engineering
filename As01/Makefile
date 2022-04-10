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

APPS = matmul_basic matmul_seq matmul_cuda transpose

.PHONY: all

all: $(APPS)

matmul_cuda: matmul_cuda.cu mmio.cpp MatrixParser/matrixparser.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

matmul_seq: matmul_seq.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

matmul_basic: matmul.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

transpose: transpose.cpp mmio.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS) *.o
