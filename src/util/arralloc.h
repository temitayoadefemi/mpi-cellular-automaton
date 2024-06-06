#ifndef __ARRALLOC_H__
#define __ARRALLOC_H__

#include <stddef.h>
#include <stdarg.h>

// Function Prototype for arralloc
// This function dynamically allocates a multi-dimensional array.
// Parameters:
//    size - size of each element in the array
//    ndim - number of dimensions
//    ...  - variable list of arguments specifying the size of each dimension
void *arralloc(size_t size, int ndim, ...);
void subarray(size_t align_size, size_t size, int ndim, int prdim, void ***pp, void **qq, int *dimp, int index);

#endif	// __ARRALLOC_H__
