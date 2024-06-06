#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "structs.h"
#include "args.h"
#include "calib.h"
#include "mplib.h"
#include "mem.h"
#include "wraplib.h"

int main(int argc, char *argv[]) {

    // Define and initialize the master structure
    master_str master;

    // Set up communication channels (MPI)
    setup_comm(&master);

    // Read and process command-line arguments
    if (read_args(&master, argc, argv) == FAILED) {

        mpstop(); // Stop the MPI environment

        return 0;
    }

    // Compute dimensions for the simulation, exit if unsuccessful
    if (compute_dimensions(&master) == FAILED) {

        mpstop(); // Stop the MPI environment

        return 0;
    }

    // Create arrays for cell data and their neighbors
    int **cell_grid = create_cell_array(&master);
    int **neighbor_grid = create_neighbours_array(&master);
    int **global_cell_grid = create_global_array(&master);
    int **reduction_cell_grid = create_reduction_array(&master);
    int **local_cell_grid = create_local_cell_array(&master);

    // Initialize and distribute the workload across the processors
    initialise_and_distribute(&master, cell_grid, global_cell_grid, local_cell_grid);

    // Process the cells based on the current simulation parameters
    process(&master, cell_grid, neighbor_grid);

    // Gather data from all nodes and write to output
    gather_write_data(&master, local_cell_grid, reduction_cell_grid, global_cell_grid, cell_grid);

    // Clean up resources and stop communication
    clean_buffers_stop_comm(&master, cell_grid, neighbor_grid, global_cell_grid, local_cell_grid, reduction_cell_grid);

    // Exit the program
    return 0;
}
