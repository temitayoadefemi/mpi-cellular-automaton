#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "structs.h"


#define NDIMS 2   
#define VERTICAL 1 
#define HORIZONTAL 0

// Initialize the MPI environment and set up the communication structure.
void mpstart(comm_str *comm) { 
    MPI_Init(NULL, NULL); // Initialize MPI environment.
    comm->comm = MPI_COMM_WORLD; // Set MPI communicator to the global communicator.
    MPI_Comm_rank(comm->comm, &comm->rank); // Get the rank of the current process.
    MPI_Comm_size(comm->comm, &comm->size); // Get the total number of processes.
}
  
// Finalize the MPI environment.
void mpstop(void) { 
    MPI_Finalize(); // Clean up all MPI state.
}

// Setup a Cartesian topology for the MPI processes.
void setup_cartesian_topology(comm_str *comm, cart_str *cart) {
    // Initialize grid dimensions and periodicity.
    cart->dims[0] = 0; // Let MPI decide.
    cart->dims[1] = 0; // Let MPI decide.
    cart->period[0] = 1; // Periodic in the first dimension.
    cart->period[1] = 0; // Non-periodic in the second dimension.
    cart->reorder = 0; // Disable reordering of processes within the grid.
    
    // Create a Cartesian grid communicator based on dimensions and periodicity.
    MPI_Dims_create(comm->size, NDIMS, cart->dims); // Automatically set the dimensions.
    MPI_Cart_create(comm->comm, NDIMS, cart->dims, cart->period, cart->reorder, &cart->comm2d); // Create the Cartesian topology.
    MPI_Cart_coords(cart->comm2d, comm->rank, NDIMS, cart->coords); // Get Cartesian coordinates of the current process.

    // Determine neighboring processes in the grid.
    MPI_Cart_shift(cart->comm2d, HORIZONTAL, VERTICAL, &cart->up.val, &cart->down.val); // Neighbors in the first dimension.
    MPI_Cart_shift(cart->comm2d, VERTICAL, VERTICAL, &cart->left.val, &cart->right.val); // Neighbors in the second dimension.
}

// Reduce the local count of cells to a global count using MPI_Reduce.
void mpi_reduce_localncell(cart_str cart, int local_live_cells, int *total_live_cells) {
    // Aggregate local cell counts across all processes.
    MPI_Reduce(&local_live_cells, total_live_cells, 1, MPI_INT, MPI_SUM, 0, cart.comm2d);
}

// Reduce local arrays to a global array using MPI_Reduce.
void mpi_reduce_allcell(cart_str cart, int **reduction_cell_grid, int **global_cell_grid, int size) {
    // Sum 2D arrays from all processes into a single global 2D array on the root process.
    MPI_Reduce(&reduction_cell_grid[0][0], &global_cell_grid[0][0], size*size, MPI_INT, MPI_SUM, 0, cart.comm2d);
}

// Broadcast the global array of cells to all processes.
void mpbcast(cart_str cart, int **global_cell_grid, int size) {
    // Distribute the global cell array from the root process to all other processes in the communicator.
    MPI_Bcast(&global_cell_grid[0][0], size*size, MPI_INT, 0, cart.comm2d);
}

// Initialize MPI data types for row and column transfers.
void initialize_mpi_types(MPI_Datatype *column_type, MPI_Datatype *row_type, master_str *master) {
    // Create a vector type for transferring columns.
    MPI_Type_vector(master->dimensions.rows, 1, master->dimensions.cols + 2, MPI_INT, column_type);
    MPI_Type_commit(column_type); // Commit the type to use it for MPI operations.

    // Create a contiguous type for transferring rows.
    MPI_Type_contiguous(master->dimensions.cols, MPI_INT, row_type);
    MPI_Type_commit(row_type); // Commit the type to use it for MPI operations.
}

// Initialize a buffer for MPI buffered send operations.
void initialize_mpi_buffer(void **buffer, int *bsize, master_str *master) {
    // Calculate the required buffer size.
    *bsize = (master->dimensions.rows + master->dimensions.cols) * sizeof(int) + MPI_BSEND_OVERHEAD;
    
    // Allocate the buffer.
    *buffer = malloc(*bsize);
    if (*buffer == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort MPI execution if memory allocation fails.
    }

    // Attach the buffer for use in buffered send operations.
    MPI_Buffer_attach(*buffer, *bsize);
}

// Send halo cells to neighboring processes.
void send_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, MPI_Request reqs[], master_str *master) {

    MPI_Isend(&cell_grid[master->dimensions.rows][1], 1, row_type, cart.down.val, 1, cart.comm2d, &reqs[0]); // Send bottom row.
    MPI_Isend(&cell_grid[1][1], 1, row_type, cart.up.val, 2, cart.comm2d, &reqs[2]); // Send top row.
    MPI_Isend(&cell_grid[1][master->dimensions.cols], 1, column_type, cart.right.val, 3, cart.comm2d, &reqs[4]); // Send right column.
    MPI_Isend(&cell_grid[1][1], 1, column_type, cart.left.val, 4, cart.comm2d, &reqs[6]); // Send left column.
}

// Receive halo cells from neighboring processes.
void receive_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, MPI_Request reqs[], master_str *master) {

    MPI_Irecv(&cell_grid[0][1], 1, row_type, cart.up.val, 1, cart.comm2d, &reqs[1]); // Receive top row.
    MPI_Irecv(&cell_grid[master->dimensions.rows+1][1], 1, row_type, cart.down.val, 2, cart.comm2d, &reqs[3]); // Receive bottom row.
    MPI_Irecv(&cell_grid[1][0], 1, column_type, cart.left.val, 3, cart.comm2d, &reqs[5]); // Receive left column.
    MPI_Irecv(&cell_grid[1][master->dimensions.cols+1], 1, column_type, cart.right.val, 4, cart.comm2d, &reqs[7]); // Receive right column.
}

// Coordinate the exchange of halo cells around the grid.
void exchange_halo_cells(int **cell_grid, MPI_Datatype row_type, MPI_Datatype column_type, cart_str cart, master_str *master) {
    MPI_Status status[8];
    MPI_Request reqs[8];

    // Initiate asynchronous sends and receives.
    send_halo_cells(cell_grid, row_type, column_type, cart, reqs, master);
    receive_halo_cells(cell_grid, row_type, column_type, cart, reqs, master);

    // Wait for all communication operations to complete.
    MPI_Waitall(8, reqs, status);
}

double mpgsum(cart_str cart, double *local_sum)
{
  double global_sum;
 
  MPI_Allreduce(local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, cart.comm2d);

  return global_sum;
} 

double gettime(void)
{ 
  return MPI_Wtime(); 
} 
