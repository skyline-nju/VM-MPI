#pragma once
#include "config.h"
#include "vect.h"
#include "comn.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

inline Vec_2<int> get_proc_rank(const Vec_2<int>& proc_size, MPI_Comm group_comm) {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  return Vec_2<int>(my_rank % proc_size.x, my_rank / proc_size.x);
#else
  return Vec_2<int>(0, 0);
#endif
}
 
/*
 * @brief Simulation domain in 2D
 * 
 * The simulation domain is usually a rectangle box, with length l_(lx, ly) and
 * origin o_(o_x, o_y). If using MPI, the domain with size gl_l_ is divided into
 * subdomains with size l_ and origin o_.
 *
 * ! The subdomains initialized by the constructor function is of the same size.
 * Use function adjust() to specify the size and origin of each subdomain.
 */
class Domain_2 {
public:
  Domain_2(const Vec_2<double>& gl_l,
           const Vec_2<int>& proc_size,
           MPI_Comm group_comm)
           : comm_(group_comm),
             proc_size_(proc_size),
             proc_rank_(get_proc_rank(proc_size, group_comm)),
             gl_l_(gl_l),
             l_(gl_l.x / proc_size.x, gl_l.y / proc_size.y),      
             o_(l_* proc_rank_) {}

  const Vec_2<double>& gl_l() const { return gl_l_; }
  const Vec_2<double>& l() const { return l_; }
  const Vec_2<double> origin() const { return o_; }
  const Vec_2<int>& proc_size() const { return proc_size_; }
  const Vec_2<int>& proc_rank() const { return proc_rank_; }
  MPI_Comm comm() const { return comm_; }

#ifdef USE_MPI
  template <typename T>
  void find_neighbor(T(*neighbor)[2]) const;
#endif

  void adjust(const Vec_2<double>& l0, const Vec_2<double>& o) { l_ = l0; o_ = o; }

  template <typename TPar>
  bool contain_particle(const TPar& p) const;

protected:
  MPI_Comm comm_;

  Vec_2<int> proc_size_;
  Vec_2<int> proc_rank_;
  Vec_2<double> gl_l_;
  Vec_2<double> l_;
  Vec_2<double> o_;
};

#ifdef USE_MPI
template <typename T>
void Domain_2::find_neighbor(T (*neighbor)[2]) const {
  const int nx = proc_size_.x;
  for (int dim = 0; dim < 2; dim++) {
    if (proc_size_[dim] > 1) {
      Vec_2<int> prev(proc_rank_);
      Vec_2<int> next(proc_rank_);
      prev[dim] -= 1;
      next[dim] += 1;
      if (prev[dim] < 0)
        prev[dim] += proc_size_[dim];
      if (next[dim] >= proc_size_[dim])
        next[dim] -= proc_size_[dim];
      neighbor[dim][0] = prev.x + prev.y * nx;
      neighbor[dim][1] = next.x + next.y * nx;
    } else {
      neighbor[dim][0] = neighbor[dim][1] = MPI_PROC_NULL;
    }
  }
}
#endif

template <typename TPar>
bool Domain_2::contain_particle(const TPar& p) const {
  const Vec_2<double> end = o_ + l_;
  return p.pos.x >= o_.x && p.pos.x < end.x&& p.pos.y >= o_.y && p.pos.y < end.y;
}

/**
 * @brief Domain with periodic boundary condition in two dimension.
 * 
 * Use flag_PBC_ to indicate whether the periodic bondary condition is applied
 * along x and y direction.
 */
class PeriodicDomain_2 : public Domain_2 {
public:
  PeriodicDomain_2(const Vec_2<double>& gl_l, 
                   const Vec_2<int>& proc_size, 
                   MPI_Comm group_comm) 
                   :  Domain_2(gl_l, proc_size, group_comm), 
                      gl_half_l_(0.5 * gl_l_), origin_(origin()),
                      flag_PBC_(proc_size_.x == 1, proc_size_.y == 1) {}

  void tangle(Vec_2<double>& pos) const;
  void untangle(Vec_2<double>& v) const;

protected:
  Vec_2<double> gl_half_l_;
  Vec_2<double> origin_;
  Vec_2<bool> flag_PBC_;
};

inline void PeriodicDomain_2::tangle(Vec_2<double>& pos) const {
  if (flag_PBC_.x) {
    tangle_1(pos.x, gl_l_.x);
  }
  if (flag_PBC_.y) {
    tangle_1(pos.y, gl_l_.y);
  }
}

inline void PeriodicDomain_2::untangle(Vec_2<double>& v) const {
  if (flag_PBC_.x) {
    untangle_1(v.x, gl_l_.x, gl_half_l_.x);
  }
  if (flag_PBC_.y) {
    untangle_1(v.y, gl_l_.y, gl_half_l_.y);
  }
}

/**
 * @brief The simulation domain is divided into grids.
 * 
 * The mesh size of the grids should be equal to or slightly larger than the
 * cutoff radius of interactions between particles. The grids described by
 * this class are real grids, while the grids in celllist class may contain
 * ghost grids.
 * 
 */
class Grid_2 {
public:
  template <typename TDomain>
  Grid_2(TDomain& dm, double r_cut);

  const Vec_2<int>& n() const { return n_; }
  const Vec_2<int>& gl_n() const { return gl_n_; }
  const Vec_2<int>& origin() const { return origin_; }
  const Vec_2<double>& lc() const { return lc_; }
  const Vec_2<double> one_over_lc() const { return Vec_2<double>(1. / lc_.x, 1. / lc_.y); }
private:
  Vec_2<int> n_;
  Vec_2<int> gl_n_;
  Vec_2<int> origin_;
  Vec_2<double> lc_; // length of each cell
};


template <typename TDomain>
Grid_2::Grid_2(TDomain& dm, double r_cut) {
  gl_n_.x = int(dm.gl_l().x / r_cut);
  gl_n_.y = int(dm.gl_l().y / r_cut);
  n_ = gl_n_;
  origin_ = Vec_2<int>(0, 0);

  lc_.x = dm.gl_l().x / gl_n_.x;
  lc_.y = dm.gl_l().y / gl_n_.y;

  std::cout << "grid size" << gl_n_ << std::endl;

#ifdef USE_MPI
  Vec_2<int> n_per_proc(gl_n_.x / dm.proc_size().x, gl_n_.y / dm.proc_size().y);
  origin_ = n_per_proc * dm.proc_rank();

  if (dm.proc_rank().x == dm.proc_size().x - 1) {
    n_.x = gl_n_.x - dm.proc_rank().x * n_per_proc.x;
  } else {
    n_.x = n_per_proc.x;
  }
  if (dm.proc_rank().y == dm.proc_size().y - 1) {
    n_.y = gl_n_.y - dm.proc_rank().y * n_per_proc.y;
  } else {
    n_.y = n_per_proc.y;
  }

  //! adjust the length and origin of subdomains
  dm.adjust(lc_ * n_, lc_ * origin_);

  int* nx = new int[dm.proc_size().x * dm.proc_size().y];
  int* ny = new int[dm.proc_size().x * dm.proc_size().y];
  MPI_Gather(&n_.x, 1, MPI_INT, nx, 1, MPI_INT, 0, dm.comm());
  MPI_Gather(&n_.y, 1, MPI_INT, ny, 1, MPI_INT, 0, dm.comm());

  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);
  if (my_rank == 0) {
    std::cout << "domain size: " << dm.gl_l().x << " * " << dm.gl_l().y << std::endl;
    std::cout << "proc size: " << dm.proc_size().x << " * " << dm.proc_size().y << std::endl;
    std::cout << "grid size in x direction:";
    for (int i = 0; i < dm.proc_size().x; i++) {
      std::cout << " " << nx[i];
    }
    std::cout << "\ngrid size in y direction:";
    for (int i = 0; i < dm.proc_size().y; i++) {
      std::cout << " " << ny[i * dm.proc_size().x];
    }
    std::cout << std::endl;
  }
  delete[] nx;
  delete[] ny;
#endif
}

/**
 * @brief Decompose the simulation domain into several subdomains.
 * 
 * @tparam T       default is double
 * @param gl_l     length of the whole domain
 * @return const Vec_2<int>  sizes of subdomain in x and y direction
 */
template <typename T>
const Vec_2<int> decompose_domain(const Vec_2<T>& gl_l, MPI_Comm group_comm) {
  Vec_2<int> proc_size{};
#ifdef USE_MPI
  int tot_proc, my_proc;

  MPI_Comm_size(group_comm, &tot_proc);
  MPI_Comm_rank(group_comm, &my_proc);
  if (gl_l.x == gl_l.y) {
    int nx = int(std::sqrt(tot_proc));
    while (tot_proc % nx != 0) {
      nx--;
    }
    proc_size = Vec_2<int>(nx, tot_proc / nx);
  } else if (gl_l.x > gl_l.y) {
    const double ratio = gl_l.x / gl_l.y;
    if (my_proc == 0) {
      std::cout << "ratio = " << ratio << std::endl;
    }
    int ny = 1;
    while (ny <= std::sqrt(tot_proc)) {
      if (ratio * ny * ny >= tot_proc && tot_proc % ny == 0) {
        proc_size = Vec_2<int>(tot_proc / ny, ny);
        break;
      }
      ny++;
    }
  } else {
    const double ratio = gl_l.y / gl_l.x;
    int nx = 1;
    while (nx <= std::sqrt(tot_proc)) {
      if (ratio * nx * nx >= tot_proc && tot_proc % nx == 0) {
        proc_size = Vec_2<int>(nx, tot_proc / nx);
        break;
      }
      nx++;
    }
  }
  if (my_proc == 0) {
    std::cout << "decompose the domain (" << gl_l.x << ", " << gl_l.y
              << ") into " << proc_size.x << " by "
              << proc_size.y << " subdomains, each of which has size "
              << gl_l.x / proc_size.x << " by " << gl_l.y / proc_size.y
              << std::endl;
  }
#else
  proc_size.x = proc_size.y = 1;
#endif
  return proc_size;
}