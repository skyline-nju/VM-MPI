#include "domain3D.h"
#include "cellList3D.h"

Vec_3<int> Domain_3::partitioning() {
  Vec_3<int> domain_size{};
#ifdef USE_MPI
  const int tot_proc = get_proc_num();
#else
  const int tot_proc = 1;
#endif

  if (tot_proc % 2 == 1 && tot_proc > 4) {
    std::cout << "tot_rank = " << tot_proc << std::endl;
    exit(1);
  }
  if (tot_proc <= 4) {
    domain_size = Vec_3<int>(1, 1, tot_proc);
  } else {
    switch (tot_proc) {
    case 6:
      domain_size = Vec_3<int>(1, 2, 3);
      break;
    case 8:
      domain_size = Vec_3<int>(2, 2, 2);
      break;
    case 10:
      domain_size = Vec_3<int>(1, 2, 5);
      break;
    case 12:
      domain_size = Vec_3<int>(2, 2, 3);
      break;
    case 24:
      domain_size = Vec_3<int>(2, 3, 4);
      break;
    case 36:
      domain_size = Vec_3<int>(3, 3, 4);
      break;
    case 48:
      domain_size = Vec_3<int>(3, 4, 4);
      break;
    case 60:
      domain_size = Vec_3<int>(3, 4, 5);
      break;
    case 72:
      domain_size = Vec_3<int>(3, 4, 6);
      break;
    case 84:
      domain_size = Vec_3<int>(3, 4, 7);
      break;
    case 96:
      domain_size = Vec_3<int>(4, 4, 6);
      break;
    default:
      std::cout << "tot_rank = " << tot_proc << std::endl;
      exit(1);
    }
  }
  return domain_size;
}

ExtDomain_3::ExtDomain_3(const Vec_3<int> &domain_size, const Vec3d & gl_l,
                         double r_cut, int gl_np)
  :size_(domain_size), gl_par_num_(gl_np) {
  find_neighbor_proc(size_, rank_, neighbor_proc_, flag_ext_);
  divide_cell(gl_l, r_cut, gl_cell_size_, cell_len_);
  divide_domain(size_, rank_, gl_cell_size_, my_cell_size_, my_first_cell_);
  my_par_num_ = divide_par_evenly(gl_par_num_);
  l_ = cell_len_ * my_cell_size_;
  half_l_ = l_ * 0.5;
  origin_ = cell_len_ * my_first_cell_;
  end_pnt_ = origin_ + l_;

  my_proc_ = get_proc_rank();
  proc_size_ = get_proc_num();
  //std::cout << "proc: " << my_proc_ << "\t" << my_par_num_ << std::endl;
}
