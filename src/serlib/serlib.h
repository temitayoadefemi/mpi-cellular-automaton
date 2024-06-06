#ifndef SERLIB_H
#define SERLIB_H

#include "structs.h"

// Initializes the communication settings for serial processing environments
void ser_initialise_comm(master_str *master);

// Reads and processes command-line arguments, returns non-zero if there is an error
int ser_read_args(int argc, char **argv, master_str *master);

// Initializes the buffers for serial computation (definition assumed to be elsewhere)
void ser_initialise_buffers(master_str *master);

// Initializes and distributes cells for serial processing
void ser_initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid, int live_cells);

// Processes cell data in a serial manner, modifying cell states based on neighbor interactions
void ser_process(master_str *master, int **cell_grid, int **neighbor_grid);

// Gathers data from serial computation and writes it to files or other outputs
void ser_gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid);

// Cleans up buffers and stops communications, preparing for shutdown in a serial environment
void ser_clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid);

// Starts timing for performance measurement in serial computation
void ser_start_timing(master_str *master);

// Stops timing and calculates the elapsed time for performance assessment
void ser_stop_timing(master_str *master);

// Prints the average timing information for the duration of the serial computation
void ser_print_timing(master_str *master);

// Applies periodic boundary conditions to cells in a serial computation environment
void ser_periodic_boundary(int **cell_grid, master_str *master);

// Sets specific boundary conditions based on the cell position and predefined boundaries
void ser_boundary_conditions(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master);

#endif // SERLIB_H
