#pragma once
#include "run2D.h"
#ifdef RANDOM_OBSTACLE
#ifdef CONTINUE_DYNAMIC

void run_RO(int gl_par_num, const Vec_2<double>& gl_l,
            double eta, double rho_s, double eps,
            unsigned long long seed, unsigned long long seed2,
            int n_step, std::string& ini_mode,
            MPI_Comm group_comm, MPI_Comm root_comm);

#endif

#endif