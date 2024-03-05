To run the metropolis algorithm using a step size of delta = 2 while varying beta:
```
cd build
cmake ..
make
./scripts/main
```

Map:

The header files are located in the `include` directory, source files are located in the `src` directory, and the implementations are located in the `scripts` directory. `newton.cxx` is supposed to implement the Newton-Raphson method for finding the value of the dimensionless `a` variable. The wavefunction and energy functional for the H2 molecule are implemented in the `src/wavefunction.cxx` file. I implemented the metropolis algorithm as a class within the `src/Metropolis.cxx` file, and the way I'm trying to use it is shown in `scripts/main.cxx` where I am varying beta while keeping all other variables constant and computing the energy.
