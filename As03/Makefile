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

APPS = matmul_cuda csc_spmv csr_spmv csr_cuda # spmv coo_spmv 

.PHONY: all

all: $(APPS)

# spmv: spmv.c mmio.c spmv_kernel.c csr_spmv_kernel.c csc_spmv_kernel.c
# 	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -o $@ $^ $(LFLAGS)

csr_spmv: spmv.c mmio.c csr_spmv_kernel.c
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -D CSR -o $@ $^ $(LFLAGS)

csc_spmv: spmv.c mmio.c csc_spmv_kernel.c
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -D CSC -o $@ $^ $(LFLAGS)

# coo_spmv: spmv.c mmio.c coo_spmv_kernel.c
# 	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -D COO -o $@ $^ $(LFLAGS)

matmul_cuda: matmul_cuda.cu mmio.cpp MatrixParser/matrixparser.cpp
	nvcc -O3 -Xptxas -O3 -o $@ $^

csr_cuda: csr.cu 
	nvcc -O3 -Xptxas -O3 -o $@ $^

.PHONY: clean
clean:
	rm -f $(APPS)
