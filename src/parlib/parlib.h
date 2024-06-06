#ifndef PARLIB_H
#define PARLIB_H

#include "structs.h"

// Initializes the MPI communication settings for parallel processing
void par_initialise_comm(master_str *master);

// Initializes the buffers for parallel computation
void par_initialise_buffers(master_str *master);

// Initializes and distributes cells for parallel processing across multiple processors
void par_initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid, int live_cells);

// Processes cell data in parallel, modifying cell states based on neighbor interactions
void par_process(master_str *master, int **cell_grid, int **neighbor_grid);

// Gathers data from parallel computation nodes and writes it to files or other outputs
void par_gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid);

// Cleans up buffers and stops communications in a parallel environment, preparing for shutdown
void par_clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid);

// Starts timing for performance measurement in parallel computation
void par_start_timing(master_str *master);

// Stops timing and calculates the elapsed time for performance assessment in parallel computation
void par_stop_timing(master_str *master);

// Prints the average timing information for the duration of the parallel computation
void par_print_timing(master_str *master);

#endif // PARLIB_H
