#ifndef MPILIB_H
#define MPILIB_H

#include <mpi.h>
#include "structs.h"  // Include necessary structures like cart_str, comm_str, etc.

// Starts the MPI environment
void mpstart(comm_str *comm);

// Stops the MPI environment, finalizing all MPI communication
void mpstop(void);

// Sets up a Cartesian topology based on the communication structure
void setup_cartesian_topology(comm_str *comm, cart_str *cart);

// Reduces local cell counts to a global count across all processes
void mpi_reduce_localncell(cart_str cart, int local_live_cells, int *total_live_cells);

// Broadcasts data to all processes in the MPI topology
void mpbcast(cart_str cart, int **global_cell_grid, int size);

// Initializes MPI data types for row and column communications
void initialize_mpi_types(MPI_Datatype *column_type, MPI_Datatype *row_type, master_str *master);

// Allocates and attaches an MPI buffer for optimized communication
void initialize_mpi_buffer(void **buffer, int *bsize, master_str *master);

// Reduces data from all cells across processes, used in gathering operations
void mpi_reduce_allcell(cart_str cart, int **reduction_cell_grid, int **global_cell_grid, int size);

// Sends boundary cell data to adjacent processes
void send_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, MPI_Request reqs[], master_str *master);

// Receives boundary cell data from adjacent processes
void receive_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, MPI_Request reqs[], master_str *master);

// Coordinates the exchange of boundary cells between adjacent processes
void exchange_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, master_str *master);

// Computes the global sum of a variable across all processes in the MPI topology
double mpgsum(cart_str cart, double *local_sum);

// Returns the current time, useful for performance measurement
double gettime(void);

#endif // MPILIB_H
