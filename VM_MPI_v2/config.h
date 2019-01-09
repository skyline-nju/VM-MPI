#pragma once

#define OUTPUT_ON
#ifdef OUTPUT_ON
#define NP_PER_NODE 12
#endif

#define SCALAR_NOISE
#ifndef SCALAR_NOISE
#define COUNT_NEIGHBOR
#endif

#define DENSITY_NOISE
#ifndef COUNT_NEIGHBOR
#define COUNT_NEIGHBOR
#endif
#define COUNT_NEIGHBOR