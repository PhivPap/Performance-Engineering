.SECONDEXPANSION:
CXX_FLAGS = -std=c++17 -O3 -g

DEF_DEPS = InputOutput/IOBody.cpp InputOutput/Parser.cpp InputOutput/Writer.cpp InputOutput/ArgParser.cpp
APPS = naive bh_naive bh_omp bh_omp_tree

all: $(APPS)

F0 = AllVersions/Naive
naive: $(F0)/naive.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

F1 = AllVersions/BH_Naive
bh_naive: $(F1)/bh_naive.cpp $(F1)/util.cpp $(F1)/quad.cpp $(F1)/body.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

F2 = AllVersions/BH_OMP
bh_omp: $(F2)/bh_omp.cpp $(F2)/util.cpp $(F2)/quad.cpp $(F2)/body.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -fopenmp -o $@ $^

F3 = AllVersions/BH_OMP_Tree
bh_omp_tree: $(F3)/bh_omp.cpp $(F3)/util.cpp $(F3)/quad.cpp $(F3)/body.cpp $(DEF_DEPS)
	$(CXX) $(CXX_FLAGS) -fopenmp -o $@ $^

clean:
	rm -f $(APPS) *.o
