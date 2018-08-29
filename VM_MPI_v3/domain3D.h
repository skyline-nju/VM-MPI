#pragma once
#define USE_MPI

#include <vector>
#include "vect.h"
#include "comn.h"

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
                    Vec3d &origin, Vec3b &flag_c);
  explicit Domain_3(const Vec3d &gl_l, const Vec3i &gl_cells_size,
                    Vec3i &cells_size, Vec3d &origin, Vec3b &flag_c);
  template <typename TNode>
  explicit Domain_3(const Vec3d &gl_l, CellListNode_3<TNode> **cl, double r_cut = 1.);
  
  void find_neighbor(Vec3i &rank, Vec3b &flag_comm, int neighbor[3][2]) const;

  void set_l(const Vec3i &gl_cells_size, Vec3i &cells_size,
             Vec3d &l, Vec3d &origin) const;

  void set_max_buf_size(int gl_par_num, double amplification);

  const Vec3d & l() const { return l_; }
  const Vec3d & gl_l() const { return gl_l_; }
  const Vec3d & gl_half_l() const { return gl_half_l_; }
  const Vec3d & origin() const { return origin_; }
  const Vec3b & flag_comm() const { return flag_comm_; }
  const Vec3i & domain_sizes() const { return gl_size_; }
  const Vec3i & domain_rank() const { return rank_; }
  const Vec3i & gl_cells_size() const { return gl_cells_size_; }
  const Vec3i & cells_size() const { return cells_size_; }
  const int & max_buf_size() const { return max_buf_size_; }
  
  static Vec_3<int> partition(const Vec3d &l, int n_proc);
  static Vec_3<int> partition(const Vec3d &l);

  int neighbor[3][2]{};
  Vec_3<block_t> inner_shell[2]{};
  Vec_3<block_t> outer_shell[2]{};
protected:
  Vec3d l_;
  Vec3d origin_{};
  Vec3d gl_l_;
  Vec3d gl_half_l_;
  
  Vec3i gl_size_{1, 1, 1}; 
  Vec3i rank_{};
  Vec3b flag_comm_{};

  Vec3i gl_cells_size_{};
  Vec3i cells_size_{};

  int max_buf_size_ = 0;
};

template <typename TNode>
Domain_3::Domain_3(const Vec3d& gl_l, CellListNode_3<TNode>** cl, double r_cut)
  : gl_l_(gl_l), gl_half_l_(gl_l * 0.5), gl_size_(partition(gl_l_)) {
  Vec_3<double> cell_len{};
  CellListNode_3<TNode>::partition(gl_l_, r_cut, gl_cells_size_, cell_len);
  find_neighbor(rank_, flag_comm_, neighbor);
  set_l(gl_cells_size_, cells_size_, l_, origin_);
  set_comm_block(cells_size_, flag_comm_, inner_shell, outer_shell);
  *cl = new CellListNode_3<TNode>(cells_size_, cell_len, gl_l_, origin_, flag_comm_);
}
