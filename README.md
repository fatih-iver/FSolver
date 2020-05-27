# FSolver

FSolver is a program that utilizes MPI to solve the systems of equations that arise from the discretization of Poisson’s equation on a cubic domain.

## Installation

There are several implementations of the MPI standard. This program is developed by using MS-MPI (v10.1.2) but also tested against Open MPI (1.10.2). 

To install Open MPI for Linux based systems:
```bash
sudo apt update
sudo apt install build-essential
sudo apt install openmpi-bin
```

## Usage


To use it on Linux based systems:

```bash
mpic++ -std=c++11 FSolver.cpp

#mpirun -np <PROCESSOR_NUMBER> a.out <n>

# PROCESSOR_NUMBER must be the cube of a positive whole number. 
# Possible values for PROCESSOR_NUMBER: 1, 8, 27, 64 ...

# n + 1 must be divisible by the cube root of the PROCESSOR_NUMBER without any remainder.

mpirun -np 8 a.out 5
```

## Conditions



## Performance Results

Test CPU: Intel(R) Core(TM) İ7-7500U CPU @ 2.70 GHZ

This processor has 2 physical and 4 logical cores.

```bash
mpic++ -std=c++11 FSolver.cpp

# No Speed Up (CPU has only 2 physical cores)
mpirun -np 1 a.out 25 # 6.1544 sn
mpirun -np 8 a.out 25 # 18.6954 sn

# No Speed Up (CPU has only 2 physical cores)
mpirun -np 1 a.out 31 # 16.0943 sn
mpirun -np 8 a.out 31 # 40.1093 sn

# No Gain (CPU has only 2 physical cores)
mpirun -np 1 a.out 47 # 132.783 sn
mpirun -np 8 a.out 47 # 130.29 sn

# Speed Up Ratio: 1.6 (CPU has only 2 physical cores)
mpirun -np 1 a.out 63 # 599.122 sn
mpirun -np 8 a.out 63 # 377.707 sn

# Speed Up Ratio: 2.4 (CPU has only 2 physical cores)
mpirun -np 1 a.out 81 # 2589.78 sn
mpirun -np 8 a.out 81 # 1100.24 sn

# Speed Up Ratio: 1.9 (CPU has only 2 physical cores)
mpirun -np 1 a.out 99 #  6407.75 sn
mpirun -np 8 a.out 99 #  3402.22 sn
```

## Justification

As you can see, there is no performance gain up to a point. The parallel program runs much slower for smaller n. When n gets bigger, performance gains can be observed. 

Be careful that, test CPU has only 2 physical cores. This is the main bottleneck and the explanation for the observed speedup values. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
