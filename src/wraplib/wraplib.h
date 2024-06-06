#ifndef __WRAPLIB_H__
#define __WRAPLIB_H__

#include "structs.h"

// Sets up the communication channels for distributed or parallel execution
void setup_comm(master_str *master);

// Parses command-line arguments and configures the master structure
status read_args(master_str *master, int argc, char **argv);

// Initializes and distributes the computational workload among available resources
void initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid);

// Executes the main processing logic based on the computation model (serial or parallel)
void process(master_str *master, int **cell_grid, int **neighbor_grid);

// Cleans up and deallocates memory buffers, stops communication channels
void clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid);

// Starts timing for performance measurement, usually used for benchmarking
void start_timing(master_str *master);

// Stops the timing, captures the duration for performance assessment
void stop_timing(master_str *master);

// Gathers data from distributed systems and prepares it for output or storage
void gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid);

#endif // __WRAPLIB_H__
