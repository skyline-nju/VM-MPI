#pragma once

#define USE_MPI

#define POLAR_ALIGN
#define SCALAR_NOISE

// #define RANDOM_TORQUE
//#define RANDOM_FIELD
#define RANDOM_OBSTACLE
#define SCATTER_NORMED

#ifndef USE_MPI
#define MPI_PROC_NULL -1
#endif

//#define INI_ORDERED
