# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -g -O0 -march=native


APPS = sync-overhead no-sync-overhead

.PHONY: all

all: $(APPS)

sync-overhead: sum.cpp 
	$(CXX) -fopenmp $(CXX_FLAGS) -D OVERHEAD -o $@ $^ $(LFLAGS)

no-sync-overhead: sum.cpp 
	$(CXX) -fopenmp $(CXX_FLAGS) -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS) *.o
