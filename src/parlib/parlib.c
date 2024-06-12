#include <string.h>
#include "calib.h"
#include "arralloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parlib.h"
#include "mem.h"
#include "misc.h"
#include "mplib.h"

#define FIRSTPERIODICBOUNDARYDIVISOR 8
#define SECONDPERIODICBOUNDARYDIVISOR 7
#define OFFSET 1

// Initializes the MPI communication and sets up the Cartesian topology for parallel computation
void par_initialise_comm(master_str *master) {
    mpstart(&master->comm);
    setup_cartesian_topology(&master->comm, &master->cart);
}

// Initializes and distributes data structures across processes for parallel computation
void par_initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid, int live_cells) {
    if (master->comm.rank == 0) {
        printf("automaton: running on %d process(es)\n", master->comm.size);
        printf("automaton: L = %d, rho = %f, seed = %d, maxstep = %d\n",
               master->params.landscape, master->params.rho, master->params.seed, master->params.maxstep);
        rinit(master->params.seed);
        initialize_cells(master->params.landscape, global_cell_grid, master, &live_cells);
    }
    mpbcast(master->cart, global_cell_grid, master->params.landscape);
    distribute_cells(local_cell_grid, global_cell_grid, master);
    copy_data_to_cell_grid(cell_grid, local_cell_grid, master);
    zero_top_bottom_halos(cell_grid, master);
    zero_left_right_halos(cell_grid, master);
}

// Processes the cell data in parallel, managing data exchange and computation across processes
void par_process(master_str *master, int **cell_grid, int **neighbor_grid) {
    MPI_Datatype column_type, row_type;
    void *buffer;
    int bsize;

    initialize_mpi_types(&column_type, &row_type, master);
    initialize_mpi_buffer(&buffer, &bsize, master);

    int local_live_cells, total_live_cells;
    par_start_timing(master);

    for (int step = 1; step <= master->params.maxstep; step++) {
        exchange_halo_cells(cell_grid, row_type, column_type, master->cart, master);
        int periodic_boundary_start = master->params.landscape / FIRSTPERIODICBOUNDARYDIVISOR + OFFSET;
        int periodic_boundary_end = (SECONDPERIODICBOUNDARYDIVISOR * master->params.landscape) / FIRSTPERIODICBOUNDARYDIVISOR;
        adjust_boundaries(cell_grid, master->cart, periodic_boundary_start, periodic_boundary_end, master);
        calculate_neighbors(cell_grid, neighbor_grid, master);
        update_cells(cell_grid, neighbor_grid, &local_live_cells, master);
        mpi_reduce_localncell(master->cart, local_live_cells, &total_live_cells);

        if (master->comm.rank == 0 && (step % master->params.printfreq == 0)) {
            printf("automaton: number of live cells on step %d is %d\n", step, total_live_cells);
            if (should_terminate(total_live_cells, master, step)) {
                break;  // Terminate if function returns true
            }
        }
    }

    par_stop_timing(master);  // Stop timing and calculate

    if (master->comm.rank == 0) {
        par_print_timing(master);  // Print the results
    }

    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);
    MPI_Buffer_detach(&buffer, &bsize);
    free(buffer);
}

// Gathers data from all processes, combines it, and writes it to a file
void par_gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid) {
    copy_data_to_local_cell_grid(cell_grid, local_cell_grid, master);
    zerotmpcell(reduction_cell_grid, master);
    gather_cells(local_cell_grid, reduction_cell_grid, master);
    mpi_reduce_allcell(master->cart, reduction_cell_grid, global_cell_grid, master->params.landscape);

    if (master->comm.rank == 0) {
        writecelldynamic("cell.pbm", global_cell_grid, master->params.landscape);
    }
}

// Cleans up and deallocates memory, stops MPI communication to prepare for shutdown
void par_clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid) {
    deallocate_arrays(cell_grid, neighbor_grid, global_cell_grid, local_cell_grid, reduction_cell_grid);
    mpstop();
}

// Starts the timing for performance evaluation in parallel processing
void par_start_timing(master_str *master) {
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.start = gettime();
#endif
}

// Stops the timing and calculates the average processing time across all processes
void par_stop_timing(master_str *master) {
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.local = gettime() - master->time.start;
    master->time.average = mpgsum(master->cart, &master->time.local) / master->comm.size;
#endif
}

// Outputs the average timing information for the duration of the parallel processing
void par_print_timing(master_str *master) {
#ifdef TIME
    printf("Average Time for %d iterations = %f\n", master->params.maxstep, master->time.average);
#endif
}
