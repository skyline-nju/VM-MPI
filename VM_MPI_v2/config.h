#ifndef CONFIG_H
#define CONFIG_H

#define USE_MPI

#define POLAR_ALIGN
#define SCALAR_NOISE

#define RANDOM_TORQUE
//#define RANDOM_FIELD
//#define RANDOM_POTENTIAL

#ifndef USE_MPI
#define MPI_PROC_NULL -1
#endif

#define REPLICAS
//#define USE_NC

// boundary conditions:
#define REF_WALL_Y  // reflecting walls at y=0 and y=Ly

//#define REF_WALL_XY // relecting walls at both x and y directions



#endif