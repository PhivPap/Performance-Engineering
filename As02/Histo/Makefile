.SECONDEXPANSION:
CXX_FLAGS = -std=c++17 -O3 -pthread

DEPS = 
APPS = histo_per_thread histo_atomic

all: $(APPS)

$(APPS): $$@.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

clean:
	rm -f $(APPS) *.o