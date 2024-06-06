#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "args.h"
#include "calib.h"
#include "mplib.h"
#include "arralloc.h"
#include "serlib.h"
#include "mem.h"
#include "misc.h"

#define FIRSTPERIODICBOUNDARYDIVISOR 8
#define SECONDPERIODICBOUNDARYDIVISOR 7
#define OFFSET 1

// Initializes communication for serial processing
void ser_initialise_comm(master_str *master) {
    mpstart(&master->comm);
    setup_cartesian_topology(&master->comm, &master->cart);
}

// Reads arguments and returns an error if the read fails
int ser_read_args(int argc, char **argv, master_str *master) {
    if (read_parameters(master, argc, argv) == 1) {
        mpstop();  // Finalize MPI before exiting
        return 1;  // Exit the program successfully
    }
    return 0;
}

// Initializes and distributes cells across the processes
void ser_initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid, int live_cells) {
    printf("automaton: running on %d process(es)\n", master->comm.size);
    printf("automaton: L = %d, rho = %f, seed = %d, maxstep = %d\n",
           master->params.landscape, master->params.rho, master->params.seed, master->params.maxstep);
    rinit(master->params.seed);
    initialize_cells(master->params.landscape, global_cell_grid, master, &live_cells);
    copy_data_to_cell_grid(cell_grid, global_cell_grid, master);
    zero_top_bottom_halos(cell_grid, master);
    zero_left_right_halos(cell_grid, master);
}

// Processes cells, calculating and updating their states
void ser_process(master_str *master, int **cell_grid, int **neighbor_grid) {
    int live_cell_count; 
    ser_start_timing(master);
    for (int step = 1; step <= master->params.maxstep; step++) {
        // Improved variable names for clarity
        int periodic_boundary_start = master->params.landscape / FIRSTPERIODICBOUNDARYDIVISOR + OFFSET;
        int periodic_boundary_end = (SECONDPERIODICBOUNDARYDIVISOR * master->params.landscape) / FIRSTPERIODICBOUNDARYDIVISOR;
        ser_periodic_boundary(cell_grid, master);
        ser_boundary_conditions(cell_grid, master->cart, periodic_boundary_start, periodic_boundary_end, master);
        calculate_neighbors(cell_grid, neighbor_grid, master);
        update_cells(cell_grid, neighbor_grid, &live_cell_count, master);
        if (step % master->params.printfreq == 0) {
            printf("automaton: number of live cells on step %d is %d\n", step, live_cell_count);
        }
        if (should_terminate(live_cell_count, master, step)) {
            break;  // Terminate if function returns true
        }
    }
    ser_stop_timing(master);  // Stop timing and calculate
    ser_print_timing(master);  // Print the results
}



// Gathers and writes data to a file
void ser_gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid) {
    copy_data_to_local_cell_grid(cell_grid, global_cell_grid, master);
    if (master->comm.rank == 0) {
        writecelldynamic("cell.pbm", global_cell_grid, master->params.landscape);
    }
}

// Cleans up and deallocates all arrays and stops communication
void ser_clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid) {
    deallocate_arrays(cell_grid, neighbor_grid, global_cell_grid, local_cell_grid, reduction_cell_grid);
    mpstop();
}

// Enforces periodic boundary conditions on the cellular grid
void ser_periodic_boundary(int **cell_grid, master_str *master) {
    for (int j = 1; j <= master->params.landscape; j++) {
        cell_grid[0][j] = cell_grid[master->params.landscape][j];
        cell_grid[master->params.landscape+1][j] = cell_grid[1][j];
    }
}

// Applies specific boundary conditions based on position
void ser_boundary_conditions(int **cell_grid, cart_str cart, int periodic_boundary_start, int periodic_boundary_end, master_str *master) {
    for (int j = 1; j <= master->params.landscape; j++) {
        if (j < periodic_boundary_start || j > periodic_boundary_end) {
            cell_grid[0][j] = 0;
            cell_grid[master->params.landscape+1][j] = 0;
        }
    }
}

// Starts timing for performance analysis
void ser_start_timing(master_str *master) {
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.start = gettime();
#endif
}

// Stops timing and calculates average time
void ser_stop_timing(master_str *master) {
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.local = gettime() - master->time.start;
    master->time.average = mpgsum(master->cart, &master->time.local) / master->comm.size;
#endif
}

// Prints the average timing information for the simulation
void ser_print_timing(master_str *master) {
#ifdef TIME
    printf("Average Time for %d iterations = %f\n", master->params.maxstep, master->time.average);
#endif
}
