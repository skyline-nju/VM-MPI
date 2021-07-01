#pragma once
#include "run2D.h"
#ifdef DILUTE_COUPLING

void run_dilute_coupling(int gl_par_num, const Vec_2<double>& gl_l,
                         double eta, double rho_s, double eps,
                         unsigned long long seed, unsigned long long seed2,
                         int n_step, std::string& ini_mode,
                         MPI_Comm group_comm, MPI_Comm root_comm);

#endif