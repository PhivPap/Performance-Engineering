# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -O3 -Wall -mavx
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

APPS = matmul_basic matmul_opt matmul_cuda transpose #check_result

.PHONY: all

all: $(APPS)

matmul_cuda: matmul_cuda.cu mmio.cpp MatrixParser/matrixparser.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

matmul_opt: matmul_opt.cpp mmio.cpp MatrixParser/matrixparser.cpp check_result.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

matmul_basic: matmul.cpp mmio.cpp MatrixParser/matrixparser.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

transpose: transpose.cpp mmio.cpp
	$(CXX) $(CXX_FLAGS) $(GENERATE) -o $@ $^ $(LFLAGS)

# check_result: check_result.cpp mmio.cpp MatrixParser/matrixparser.cpp
# 	$(CXX) $(CXX_FLAGS) $(GENERATE) -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS) *.o
