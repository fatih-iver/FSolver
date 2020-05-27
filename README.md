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

## 6 Points Per Edge ##
mpirun -np 1 a.out 5 -> 0.0339636 sn # Base
mpirun -np 8 a.out 5 -> 0.134244 sn # SpeedUp: 0.253
mpirun -np 27 a.out 5 -> 0.454991 sn # SpeedUp: 0.075

## 30 Points Per Edge ##
mpirun -np 1 a.out 29 -> 1.51569 sn # Base
mpirun -np 8 a.out 29 -> 0.872383 sn # SpeedUp: 1.737
mpirun -np 27 a.out 29 -> 1.35544 sn # SpeedUp: 1.118

## 6O Points Per Edge ##
mpirun -np 1 a.out 59 -> 70.3674 sn # Base
mpirun -np 8 a.out 59 -> 24.7332 sn # SpeedUp: 2.845
mpirun -np 27 a.out 59 -> 26.3868 sn # SpeedUp: 2.667

## 90 Points Per Edge ##
mpirun -np 1 a.out 89 -> 560.749 sn # Base
mpirun -np 8 a.out 89 -> 187.577 sn # SpeedUp: 2.989
mpirun -np 27 a.out 89 -> 196.726 sn # SpeedUp: 2.850

```

## Justification

As you can see, there is no performance gain up to a point. The parallel program runs much slower for smaller n. When n gets bigger, performance gains can be observed. Speed up is close to 3 beacause the test CPU has 2 physical cores and 4 logical cores. Logical cores are not one to one corresponding of physical cores. In average, with hyper-threading, a machine perform %30 faster. This explains why speed up is around 3. It would be around 2 with no hyper-threading - no logical cores. Also, note that, speed up values for 8 processors and 27 processors are close. This is because they run on the same hardware even though the given processors numbers are different. However, one with the 27 processors performs a little bit worse because communication overhead increases with the number of processors.  

Be careful that, test CPU has only 2 physical cores. This is the main bottleneck and the explanation for the observed speedup values. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
