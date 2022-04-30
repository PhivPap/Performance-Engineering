.SECONDEXPANSION:
CXX_FLAGS = -std=c++17 -O3

DEPS = 
APPS = histo_atomic # histo_no_atomic ...


all: $(APPS)

$(APPS): $$@.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

clean:
	rm -f $(APPS) *.o