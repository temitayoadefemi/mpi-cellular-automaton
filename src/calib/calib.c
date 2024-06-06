#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "arralloc.h"
#include "misc.h"
#include <stdbool.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "arralloc.h"
#include "misc.h"
#include <stdbool.h>

// Compute the dimensions of the grid based on the landscape size and Cartesian grid dimensions.
int check_divisibility(int landscape, int dimension, const char* dim_name, int rank) {
    if (landscape % dimension != 0) {
        if (rank == 0) {
            printf("Error: Landscape is not perfectly divisible by %s dimension, use a number of processes that is a factor of the landscape value\n", dim_name);
        }
        return 0; // Exit if not divisible
    }
    return landscape / dimension; // Return the computed dimension
}

// Refactored function
int compute_dimensions(master_str *master) {
    if (master->params.version == par2D) {
        int LX = check_divisibility(master->params.landscape, master->cart.dims[0], "first", master->comm.rank);
        if (!LX) return FAILED;

        int LY = check_divisibility(master->params.landscape, master->cart.dims[1], "second", master->comm.rank);
        if (!LY) return FAILED;

        master->dimensions.rows = LX;
        master->dimensions.cols = LY;
    }
    else if (master->params.version == serial) {
        master->dimensions.rows = master->params.landscape;
        master->dimensions.cols = master->params.landscape;
    }

    return SUCCESS;
}

// Distribute cells from the global landscape to a smaller grid based on the process's subdomain.
void distribute_cells(int **local_cell_grid, int **global_cell_grid, master_str *master) {
    // Loop through each cell in the smaller grid to assign values from the global grid.
    for (int i = 0; i < master->dimensions.rows; i++) {
        for (int j = 0; j < master->dimensions.cols; j++) {
            int global_row = master->cart.coords[0] * master->dimensions.rows + i;
            int global_col = master->cart.coords[1] * master->dimensions.cols + j;
            local_cell_grid[i][j] = global_cell_grid[global_row][global_col];
        }
    }
}

// Set all elements of the temporary cell array to zero.
void zerotmpcell(int **reduction_cell_grid, master_str *master) {
    // Loop through each element in the temporary grid and set its value to zero.
    for (int i = 0; i < master->params.landscape; i++) {
        for (int j = 0; j < master->params.landscape; j++) {
            reduction_cell_grid[i][j] = 0;
        }
    }
}

// Copy data from the padded cell array to the smaller cell array.
void copy_data_to_local_cell_grid(int **cell_grid, int **local_cell_grid, master_str *master) {
    // Loop through each cell in the padded grid and copy it to the corresponding position in the smaller grid.
    for (int i = 1; i <= master->dimensions.rows; i++) {
        for (int j = 1; j <= master->dimensions.cols; j++) {
            local_cell_grid[i - 1][j - 1] = cell_grid[i][j];
        }
    }
}

// Gather cells from the smaller grids of each process into the global temporary grid.
void gather_cells(int **local_cell_grid, int **reduction_cell_grid, master_str *master) {
    // Loop through each cell in the smaller grid and place its value in the correct position in the global grid.
    for (int i = 0; i < master->dimensions.rows; i++) {
        for (int j = 0; j < master->dimensions.cols; j++) {
            int global_row = master->cart.coords[0] * master->dimensions.rows + i;
            int global_col = master->cart.coords[1] * master->dimensions.cols + j;
            reduction_cell_grid[global_row][global_col] = local_cell_grid[i][j];
        }
    }
}

// Initialize the cells based on a probability and count the number of live cells.
void initialize_cells(int landscape, int **global_cell_grid, master_str *master, int *live_cells) {
    int initial_live_cells = 0;  // Initialize the live cell count.
    double r;  // Variable for random probability.

    // Loop through each cell in the landscape.
    for (int i = 0; i < landscape; i++) {
        for (int j = 0; j < landscape; j++) {
            r = uni();  // Generate a random double between 0.0 and 1.0.

            // Set cell state based on probability and update live cell count.
            if (r < master->params.rho) {
                global_cell_grid[i][j] = 1;
                initial_live_cells++;
            } else {
                global_cell_grid[i][j] = 0;
            }
        }
    }

    *live_cells = initial_live_cells;  // Update the live cell count.
    // Print the density information.
    master->initialcells = initial_live_cells;
    printf("automaton: rho = %f, live cells = %d, actual density = %f\n",
           master->params.rho, initial_live_cells, ((double) initial_live_cells) / (landscape * landscape));
}

// Adjust the top boundary condition for cells, setting values based on Cartesian grid position.
void adjust_top_boundary(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master) {
    // Check if the process is at the top boundary of the Cartesian grid.
    if (cart.coords[0] == 0) {
        // Loop through cells at the top boundary.
        for (int i = 1; i < master->dimensions.cols + 1; i++) {
            int index = cart.coords[1] * master->dimensions.cols + i;
            // Adjust cells not in the periodic boundary condition range.
            if (index < periodic_boundary_start || index > periodic_boundary_end) {
                cell_grid[0][i] = 0;
            }
        }
    }
}

// Adjust the bottom boundary condition for cells, setting values based on Cartesian grid position.
void adjust_bottom_boundary(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master) {
    // Check if the process is at the bottom boundary of the Cartesian grid.
    if (cart.coords[0] == cart.dims[0] - 1) {
        // Loop through cells at the bottom boundary.
        for (int i = 1; i < master->dimensions.cols + 1; i++) {
            int index = cart.coords[1] * master->dimensions.cols + i;
            // Adjust cells not in the periodic boundary condition range.
            if (index < periodic_boundary_start || index > periodic_boundary_end) {
                cell_grid[master->dimensions.rows + 1][i] = 0;
            }
        }
    }
}

// Adjust both top and bottom boundary conditions for cells.
void adjust_boundaries(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master) {
    // Call functions to adjust top and bottom boundaries.
    adjust_top_boundary(cell_grid, cart, periodic_boundary_start, periodic_boundary_end, master);
    adjust_bottom_boundary(cell_grid, cart, periodic_boundary_start, periodic_boundary_end, master);
}

// Update each cell based on its neighbors' states and count the number of live cells.
void update_cells(int **cell_grid, int **neighbor_grid, int *local_live_cells, master_str *master) {
    *local_live_cells = 0;  // Reset the count of live cells.
    for (int i = 1; i <= master->dimensions.rows; i++) {
        for (int j = 1; j <= master->dimensions.cols; j++) {
            // If the cell has 2, 4, or 5 neighbors, it becomes or remains alive; otherwise, it dies.
            if (neighbor_grid[i][j] == 2 || neighbor_grid[i][j] == 4 || neighbor_grid[i][j] == 5) {
                cell_grid[i][j] = 1;
                (*local_live_cells)++;  // Increment live cell count.
            } else {
                cell_grid[i][j] = 0;
            }
        }
    }
}

// Calculate the number of neighbors for each cell in the grid.
void calculate_neighbors(int **cell_grid, int **neighbor_grid, master_str *master) {
    for (int i = 1; i <= master->dimensions.rows; i++) {
        for (int j = 1; j <= master->dimensions.cols; j++) {
            // Sum the states of the cell and its immediate neighbors to get the total number of active neighbors.
            neighbor_grid[i][j] = cell_grid[i][j] + cell_grid[i-1][j] + cell_grid[i+1][j] + cell_grid[i][j-1] + cell_grid[i][j+1];
        }
    }
}

// Copy data from the smaller cell array back to the main cell array after calculations.
void copy_data_to_cell_grid(int **cell_grid, int **local_cell_grid, master_str *master) {
    for (int i = 1; i <= master->dimensions.rows; i++) {
        for (int j = 1; j <= master->dimensions.cols; j++) {
            // Transfer data from local_cell_grid (local grid) back to cell_grid (global grid with ghost rows and columns).
            cell_grid[i][j] = local_cell_grid[i - 1][j - 1];
        }
    }
}

// Set halo cells to zero along the top and bottom boundaries of the grid to manage boundary conditions.
void zero_top_bottom_halos(int **cell_grid, master_str *master) {
    for (int i = 0; i <= master->dimensions.rows + 1; i++) {
        // Reset the top and bottom halo cells to zero.
        cell_grid[i][0] = 0;
        cell_grid[i][master->dimensions.cols + 1] = 0;
    }
}

// Set halo cells to zero along the left and right boundaries of the grid to manage boundary conditions.
void zero_left_right_halos(int **cell_grid, master_str *master) {
    for (int j = 0; j <= master->dimensions.cols + 1; j++) {
        // Reset the left and right halo cells to zero.
        cell_grid[0][j] = 0;
        cell_grid[master->dimensions.rows + 1][j] = 0;
    }
}

bool should_terminate(int ncell, master_str *master, int step) {
    if (ncell < 0.75 * master->initialcells || ncell > 1.33 * master->initialcells) {
        printf("Terminating early: number of live cells out of threshold range on step %d\n", step);
        return true;  // Conditions met, suggest termination
    }
    return false;  // Conditions not met, continue simulation
}
