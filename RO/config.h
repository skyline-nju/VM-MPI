#pragma once

//#define USE_MPI

#define POLAR_ALIGN
#define SCALAR_NOISE
#define CONTINUE_DYNAMIC

#define RANDOM_OBSTACLE
#define SCATTER_NORMED

//#define DILUTE_COUPLING

#ifndef USE_MPI
#define MPI_PROC_NULL -1
typedef int MPI_Comm;
#endif
