#ifndef CALIB_H
#define CALIB_H

#include "structs.h"  // Including necessary structures like cart_str, master_str, etc.
#include <stdbool.h>

// Distributes cell data from a large cell array to smaller, more manageable arrays
void distribute_cells(int **local_cell_grid, int **global_cell_grid, master_str *master);

// Zeroes out the temporary cell array
void zerotmpcell(int **reduction_cell_grid, master_str *master);

// Copies data from a general cell array to a smaller cell array
void copy_data_to_local_cell_grid(int **cell_grid, int **local_cell_grid, master_str *master);

// Gathers processed cell data from smaller arrays into a temporary array
void gather_cells(int **local_cell_grid, int **reduction_cell_grid, master_str *master);

// Initializes cell data for a given landscape, setting live cells based on parameters
void initialize_cells(int landscape, int **global_cell_grid, master_str *master, int *live_cells);

// Adjusts the top boundary of the cell grid based on predefined boundary conditions
void adjust_top_boundary(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master);

// Adjusts the bottom boundary of the cell grid based on predefined boundary conditions
void adjust_bottom_boundary(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master);

// Adjusts the top and bottom boundaries of the cell grid
void adjust_boundaries(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master);

// Updates cell states based on neighbor data
void update_cells(int **cell_grid, int **neighbor_grid, int *local_live_cells, master_str *master);

// Calculates the number of neighboring live cells for each cell in the array
void calculate_neighbors(int **cell_grid, int **neighbor_grid, master_str *master);

// Copies data from a smaller cell array back to the main cell array
void copy_data_to_cell_grid(int **cell_grid, int **local_cell_grid, master_str *master);

// Clears data in the top and bottom halo regions of the cell grid
void zero_top_bottom_halos(int **cell_grid, master_str *master);

// Clears data in the left and right halo regions of the cell grid
void zero_left_right_halos(int **cell_grid, master_str *master);

// Applies periodic boundary conditions to the cell grid
void periodic_boundary(int **cell_grid, master_str *master);

// Computes the dimensions for the cell grid based on the master settings
int compute_dimensions(master_str *master);

// Terminates the calculation if the grid exceeds or decreases past a threshold
bool should_terminate(int ncell, master_str *master, int step);

#endif // CALIB_H
