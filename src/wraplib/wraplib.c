#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "args.h"
#include "wraplib.h"
#include "serlib.h"
#include "parlib.h"
#include "structs.h"

// Initializes the communication channels based on the type of parallelization
void setup_comm(master_str *master) {
    par_initialise_comm(master);
}

// Reads command-line arguments and updates the master structure accordingly
status read_args(master_str *master, int argc, char **argv) {
    if (read_parameters(master, argc, argv) == 1) {
        return FAILED;
    }
    return SUCCESS;
}

// Initializes and distributes cells based on the execution mode (parallel or serial)
void initialise_and_distribute(master_str *master, int **cell_grid, int **global_cell_grid, int **local_cell_grid) {
    int live_cells = 0;
    if (master->params.version == par2D) {
        par_initialise_and_distribute(master, cell_grid, global_cell_grid, local_cell_grid, live_cells);
    } else if (master->params.version == serial) {
        ser_initialise_and_distribute(master, cell_grid, global_cell_grid, local_cell_grid, live_cells);
    }
}

// Processes cells based on the execution mode
void process(master_str *master, int **cell_grid, int **neighbor_grid) {
    if (master->params.version == par2D) {
        par_process(master, cell_grid, neighbor_grid);
    } else if (master->params.version == serial) {
        ser_process(master, cell_grid, neighbor_grid);
    }
}

// Starts the timing for performance measurement
void start_timing(master_str *master) {
    if (master->params.version == par2D) {
        par_start_timing(master);
    } else if (master->params.version == serial) {
        ser_start_timing(master);
    }
}

// Stops the timing for performance measurement
void stop_timing(master_str *master) {
    if (master->params.version == par2D) {
        par_stop_timing(master);
    }
    if (master->params.version == serial) {
        ser_stop_timing(master);
    }
}

// Gathers data from worker nodes and writes it to files or other outputs
void gather_write_data(master_str *master, int **local_cell_grid, int **reduction_cell_grid, int **global_cell_grid, int **cell_grid) {
    if (master->params.version == par2D) {
        par_gather_write_data(master, local_cell_grid, reduction_cell_grid, global_cell_grid, cell_grid);
    }
    if (master->params.version == serial) {
        ser_gather_write_data(master, local_cell_grid, reduction_cell_grid, global_cell_grid, cell_grid);
    }
}

// Cleans up buffers and stops communication, preparing for shutdown
void clean_buffers_stop_comm(master_str *master, int **cell_grid, int **neighbor_grid, int **global_cell_grid, int **local_cell_grid, int **reduction_cell_grid) {
    if (master->params.version == par2D) {
        par_clean_buffers_stop_comm(master, cell_grid, neighbor_grid, global_cell_grid, local_cell_grid, reduction_cell_grid);
    } else if (master->params.version == serial) {
        ser_clean_buffers_stop_comm(master, cell_grid, neighbor_grid, global_cell_grid, local_cell_grid, reduction_cell_grid);
    }
}
