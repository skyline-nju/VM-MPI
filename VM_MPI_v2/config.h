#ifndef CONFIG_H
#define CONFIG_H

#define USE_MPI

#define POLAR_ALIGN
#define SCALAR_NOISE

#define RANDOM_TORQUE
//#define RANDOM_FIELD

#ifndef USE_MPI
#define MPI_PROC_NULL -1
#endif

#define REPLICAS
//#define USE_NC
#endif