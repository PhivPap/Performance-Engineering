.SECONDEXPANSION:
CXX_FLAGS = -std=c++17 -O3 -g

DEPS = InputOutput/IOBody.cpp InputOutput/Parser.cpp InputOutput/Writer.cpp
APPS = naive # bh_naive opt0 opt1 opt2 ...


all: $(APPS)

$(APPS): AllVersions/$$@.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

clean:
	rm -f $(APPS) *.o
