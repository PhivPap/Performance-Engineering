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

APPS = csr_spmv csc_spmv

.PHONY: all

all: $(APPS)

csr_spmv: spmv.c mmio.c csr_spmv_kernel.c
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -D CSR -o $@ $^ $(LFLAGS)

csc_spmv: spmv.c mmio.c csc_spmv_kernel.c
	$(CXX) $(CXX_FLAGS) $(GENERATE) $(TIMING) -D CSC -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS)
