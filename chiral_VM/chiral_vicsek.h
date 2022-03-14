#pragma once
#include "run2D.h"

#ifdef CONTINUE_DYNAMIC

void run(int gl_par_num, const Vec_2<double>& gl_l,
            double D_0, double omega0,
            unsigned long long seed,
            int n_step, std::string& ini_mode,
            MPI_Comm group_comm, MPI_Comm root_comm);

#endif