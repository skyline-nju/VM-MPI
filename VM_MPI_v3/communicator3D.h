#pragma once

#include "vect.h"

Vec_3<int> divide_cubic_domain();

void find_neighbor_proc(const Vec_3<int> &domain_size,
                        Vec_3<int> &domain_rank,
                        int neighbor_proc[3][2],
                        Vec_3<bool> &flag_ext,
                        int test_rank = -1);

void divide_domain(const Vec_3<int>& domain_size,
                   const Vec_3<int>& domain_rank,
                   const Vec_3<int>& gl_cell_num,
                   Vec_3<int>& my_cell_num,
                   Vec_3<int>& my_first_cell);