#pragma once
#define USE_MPI
#define NP_PER_NODE 12

#define POLAR_ALIGN
#define SCALAR_NOISE

#define RANDOM_TORQUE

#define OUTPUT_ON


#ifndef USE_MPI
#define MPI_PROC_NULL -1
#endif