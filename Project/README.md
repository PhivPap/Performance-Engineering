# Peformance Engineering: 2D Barnes-Hut Approximation
University of Amsterdam - June 2022

## Authors: 
- Foivos Papapanagiotakis-Bousy
- Andre Palheiros da Silva
- Tim van Kemenade
- Dymitr Lubczy


## Usage 
### Requirements:
- g++ (GCC) v6 or higher
- GNU Make
- python3 (only for new input dataset generation)


### Installation & Compilation
```
git clone "https://github.com/PhivPap/Performance-Engineering"
cd Performance-Engineering/Project
make
```

### Input
All input files (files containg body positions, velocities and mass) can be found in *Project/Bodyfiles/in*. 

The python script which generates these files is at *Project/Bodyfiles/UniverseGenerator.py*. It does not accept any command line arguments, its configuration can be set inline.


### Execution 
The makefile will generate 4 executables: 

- *naive* refers to the O(n^2) version (prototype 0)
- *bh_naive* refers to prototype 1
- *bh_omp* refers to prototype 2
- *bh_omp_tree* refers to prototype 3

Any of the above can be run without command line arguments and the program should run in default configuration.

To customize the configuration the following flags can be set:

- ``-in`` - specifies input file
- ``-out`` - specifies output file
- ``-it`` - specifies the number of total iterations
- ``-it_len`` - specifies the number of seconds passing in between two iterations
- ``-theta`` - specifies the approximation factor theta
- ``-threads`` - specifies the number of threads used in parallel regions


Some sample executions:
```
./naive -in BodyFiles/in/in100.tsv -it 10000
./bh_naive -in BodyFiles/in/in1000.tsv -out tmp.tsv -it 200 -it_len 7200
./bh_omp -it 200 -threads 8
./bh_omp_tree -theta 1.2 -threads 16
```

**Note:** Before running the parallel versions (*bh_omp*, *bh_omp_tree*), it is recommended to pin the OpenMP threads to the cores of the processor. On linux, this can be done with:
``export OMP_PROC_BIND=true``


### Simulation Error Measurement
To measure the error produced by the Barnes-Hut approximation, the python script at *Project/Bodyfiles/ErrorCalculator.py* can be used. It receives no command line arguments. Its configuration only consists of two files; the reference output and the approximation output which can both be set inline.
