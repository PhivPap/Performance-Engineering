.SECONDEXPANSION:
CXX_FLAGS = -std=c++17 -O3

DEF_DEPS = InputOutput/IOBody.cpp InputOutput/Parser.cpp InputOutput/Writer.cpp
APPS = naive bh_naive

all: $(APPS)

F0 = AllVersions/Naive
naive: $(F0)/naive.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

F1 = AllVersions/BH_Naive
bh_naive: $(F1)/bh_naive.cpp $(F1)/util.cpp $(F1)/quad.cpp $(F1)/body.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

clean:
	rm -f $(APPS) *.o
