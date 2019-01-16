#pragma once
#define USE_MPI

//#define POLAR_ALIGN

#define OUTPUT_ON
#ifdef OUTPUT_ON
#define NP_PER_NODE 12
#endif

#define SCALAR_NOISE
#ifndef SCALAR_NOISE
#define COUNT_NEIGHBOR
#endif

//#define DISORDER_ON
#ifdef DISORDER_ON
//#define RANDOM_TORQUE 1
//#define RANDOM_FIELD 2
#define RANDOM_STRESS 3
#ifndef RANDOM_TORQUE
#define COUNT_NEIGHBOR
#endif
#endif


//#define DENSITY_NOISE
#ifdef DENSITY_NOISE
#define COUNT_NEIGHBOR
#endif

#define BIRTH_DEATH 1