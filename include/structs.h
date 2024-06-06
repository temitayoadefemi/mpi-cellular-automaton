#ifndef __STRUCTS_H__
#define __STRUCTS_H__

#include <mpi.h>
#define ndims 2 

typedef enum version_enum
{
	serial,
	par2D,

}version;

typedef struct dimensions_struct
{
	int rows;
	int cols;

} dim_str;

typedef enum {
    SUCCESS,
    FAILED
} status;



/* Contains the inforamtion about a neighbour */
typedef struct dim_comm_struct
{
	int val;

} dir_str;

/* Variables for 2D Cartesian topology */
typedef struct cartesian_struct
{
	MPI_Comm comm2d;
	int disp;
	int dims[ndims];
	int period[ndims];
	int coords[ndims];
	int reorder;
	dir_str right, left, down, up;

} cart_str;

/* Communication variables */
typedef struct communication_struct
{
	MPI_Comm comm;
	int rank, size;

} comm_str;


typedef struct time_struct
{
	double start;
	double local;
	double average;

}time_str;



typedef struct parameters
{
	  int seed;
      double rho;
	  int printfreq;
	  int landscape;
	  int maxstep;
	  double r;
	  version version;
} params_str;


typedef struct master {

    params_str params;
    comm_str comm;
    cart_str cart;
    dim_str dimensions;
    int initialcells;
    int version;
	time_str time;
} master_str;


#endif	//__STRUCTS_H__