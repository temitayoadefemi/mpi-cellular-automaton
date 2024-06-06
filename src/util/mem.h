#ifndef MEM_H
#define MEM_H

#include "structs.h"  // Including necessary structures like master_str for context

// Function declarations for memory management related to cellular automaton arrays:

// Creates a primary array for cell management in the simulation
int** create_cell_array(master_str *master);

// Creates an array for storing neighbor information of cells
int** create_neighbours_array(master_str *master);

// Creates an array for managing smaller, often sub-processed cells
int** create_local_cell_array(master_str *master);

// Creates an array that encompasses all cells in the simulation, used for comprehensive processes
int** create_global_array(master_str *master);

// Creates a temporary array used for intermediate calculations or storage
int** create_reduction_array(master_str *master);

// Deallocates all dynamic memory allocated for arrays used in the simulation
void deallocate_arrays(void *cell, void *neigh, void *allcell, void *smallcell, void *tmpcell);

// Allocates all dynamic memory allocated for arrays used in the simulation
int** allocate_2d_array(int rows, int cols);

// Handle memory allocation failure.
void handle_allocation_failure();

#endif // MEM_H




