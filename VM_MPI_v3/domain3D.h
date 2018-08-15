#pragma once

#include <vector>
#include "vect.h"
#include "comn.h"
#include "Domain3D.h"
#define USE_MPI

#ifdef USE_MPI
#include "communicator3D.h"
#endif

class Domain_3 {
public:
  typedef Vec_3<double> Vec3d;
  typedef Vec_3<int> Vec3i;
  typedef Vec_3<bool> Vec3b;

  explicit Domain_3(const Vec3d &gl_l);
  explicit Domain_3(const Vec3d &gl_l, const Vec3i &gl_size,
                    const Vec3i &gl_cells_size, Vec3i &cells_size,
                    Vec3d &origin, Vec3b &flag_comm);
  explicit Domain_3(const Vec3d &gl_l, const Vec3i &gl_cells_size,
                    Vec3i &cells_size, Vec3d &origin, Vec3b &flag_comm);
  template <typename TNode>
  explicit Domain_3(const Vec3d &gl_l, CellListNode_3<TNode> **cl, double r_cut = 1.);
  

  void find_neighbor(Vec3i &rank, Vec3b &flag_comm, int neighbor[3][2]) const;

  void set_l(const Vec3i &gl_cells_size, Vec3i &cells_size,
             Vec3d &l, Vec3d &origin) const;

  const Vec3d & gl_l() const { return gl_l_; }
  const Vec3d & gl_half_l() const { return gl_half_l_; }
  const Vec3d & origin() const { return origin_; }
  const Vec3b & flag_comm() const { return flag_comm_; }
  
  static Vec_3<int> partition(const Vec3d &l, int n_proc);
protected:
  Vec3d l_;
  Vec3d origin_{};
  Vec3d end_pnt_;

  Vec3d gl_l_;
  Vec3d gl_half_l_;
  
  Vec3i gl_size_{1, 1, 1};
  Vec3i rank_{};
  Vec3b flag_comm_{};

  int neighbor_[3][2]{};
};

template <typename TNode>
Domain_3::Domain_3(const Vec3d& gl_l, CellListNode_3<TNode>** cl, double r_cut)
  : gl_l_(gl_l), gl_half_l_(gl_l * 0.5) {
  gl_size_ = partition(gl_l_, get_proc_num());
  Vec_3<int> gl_cells_size{};
  Vec_3<int> cells_size{};
  Vec_3<double> cell_len{};
  cells_partition(gl_l_, r_cut, gl_cells_size, cell_len);
  find_neighbor(rank_, flag_comm_, neighbor_);
  set_l(gl_cells_size, cells_size, l_, origin_);
  end_pnt_ = origin_ + l_;
  *cl = new CellListNode_3<TNode>(cells_size, cell_len, gl_l_, origin_, flag_comm_);
}
