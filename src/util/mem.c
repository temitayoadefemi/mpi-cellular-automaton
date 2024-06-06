#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "arralloc.h"


#define HALO 1


// Handle memory allocation failure.
void handle_allocation_failure() {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
}

// Define a function to allocate memory for a 2D array.
int** allocate_2d_array(int rows, int cols) {
    int **array = (int**) arralloc(sizeof(int), 2, rows, cols);
    if (array == NULL) {
        handle_allocation_failure();
    }
    return array;
}


int** create_cell_array(master_str *master) {
    return allocate_2d_array(master->dimensions.rows + (HALO*2), master->dimensions.cols + (HALO*2));
}

int** create_neighbours_array(master_str *master) {
    return allocate_2d_array(master->dimensions.rows + (HALO*2), master->dimensions.cols + (HALO*2));
}

int** create_local_cell_array(master_str *master) {
    return allocate_2d_array(master->dimensions.rows, master->dimensions.cols);
}

int** create_global_array(master_str *master) {
    return allocate_2d_array(master->params.landscape, master->params.landscape);
}

int** create_reduction_array(master_str *master) {
    return allocate_2d_array(master->params.landscape, master->params.landscape);
}

void deallocate_arrays(void *cell_grid, void *neighbors_grid, void *global_cell_grid, void *local_cell_grid, void *reduction_cell_grid) {
    // Check each array pointer and free memory if it is not NULL.
    if (cell_grid) {
        free(cell_grid);
    }
    if (neighbors_grid) {
        free(neighbors_grid);
    }
    if (global_cell_grid) {
        free(global_cell_grid);
    }
    if (local_cell_grid) {
        free(local_cell_grid);
    }
    if (reduction_cell_grid) {
        free(reduction_cell_grid);
    }
}
