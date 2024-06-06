# Conway's Game of Life (Cellular Automaton)
This project implements a two-dimensional Cellular Automaton using domain decomposition and non-blocking communications, suitable for both serial and parallel execution via message-passing programming. The code implements new boundary conditions for both the serial version and the parallel version of the code. 

## What is included
- `include/`: Contains the header file called `structs.h`. This contains all the derived data structures used in the development of the code.
- `src/calib/`: Contains all the functions used to perform the cellular automaton.
- `src/mplib/`: Contains all the functions used to parallelize the code using message-passing programming.
- `src/parlib/`: Contains all the wrap functions used to generate the parallel version of the project.
- `src/serlib/`: Contains all the wrap functions used to generate the serial version of the the project.
- `src/util/`: Contains all the helper functions used in the project.
	- `args.h`: Functions that parse the command line input in the project and obtain the desired parameters and file names.
	- `arralloc.h`: Provided file that contains a function to declare an N-dimensional array avoiding the problems occuring by `malloc`.
	- `mem.h`: Contains functions that allocate and deallocate desired buffers for each implementation. Also, a function that swaps pointers to avoid copying data in each buffer.
	- `misc.h`: Contains functions that write back the data in a `.pbm` file from the buffers and also the uni and rand functions


## Options
Before compile and execute the code it is useful to know that different versions of the code can be compiled by specifying the `DEFINE` variable at the top of the `MAKEFILE`.

Available options are:

- `-DTIME`: is defined when the main loop needs to be timed.

Comment out accordingly which ones don't want to be used create a clean directory and recompile the code as it will be explained below.

Both the serial and parallel code come with the option to specify input arguments to the program through the command line. The available options are:

- `-rho`: The density of the automaton. If not specified, the default value is `0.51`.
- `-printfreq`: The frequency at which output is printed. The default frequency is `500`.
- `-landscape`: The size of the landscape to be used in the simulation. The default size is `1152`.
- `-maxstep`: The maximum number of simulation steps to be executed. The default is `10 * 1152` steps, calculated as ten times the landscape size.
- `<seed>`: The seed for the random number generator. This is a mandatory argument and must be the first argument provided.

## Usage

### Building
Select the desired `DEFINE` options from the `MAKEFILE` to compile the code with.

To run this code on cirrus, you would need to load the necessary modules:

```sh
$ module load mpt
```

To compile the serial and parallel code use:

```sh
$ make all
```

### Cleaning
To clean the project run:
```sh
$ make clean
```

### Running
To execute the serial code:
```sh

$ mpirun -n 1 `./automaton <seed> [-rho value] [-printfreq value] [-landscape value] [-maxstep value]` 

or 

$ `./automaton <seed> [-rho value] [-printfreq value] [-landscape value] [-maxstep value]` 
```

To execute the parallel code:
```sh

$ mpirun -n <int> `./automaton <seed> [-rho value] [-printfreq value] [-landscape value] [-maxstep value]` 

```