# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -g -O0 -march=native


APPS = false-sharing-on false-sharing-off

.PHONY: all

all: $(APPS)

false-sharing-on: false-sharing.c
	$(CXX) -fopenmp $(CXX_FLAGS) -D FALSESHARING -o $@ $^ $(LFLAGS)

false-sharing-off: false-sharing.c
	$(CXX) -fopenmp $(CXX_FLAGS) -o $@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -f $(APPS) *.o
