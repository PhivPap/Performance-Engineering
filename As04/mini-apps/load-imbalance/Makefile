# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -g -O0 -march=native


APPS = load-imbalance load-balance

.PHONY: all

all: $(APPS)

load-imbalance: load-balancing.cpp 
	$(CXX) -fopenmp $(CXX_FLAGS) -D IMBALANCE -o $@ $^ $(LFLAGS)

load-balance: load-balancing.cpp 
	$(CXX) -fopenmp $(CXX_FLAGS) -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS) *.o
