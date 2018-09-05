#pragma once
#define USE_MPI

#include <vector>
#include "vect.h"
#include "comn.h"

#ifdef USE_MPI
#include "communicator2D.h"
#endif

class Domain_2 {
public:
  typedef Vec_2<double> Vec2d;
  typedef Vec_2<int> Vec2i;
  typedef Vec_2<bool> Vec2b;

  explicit Domain_2(const Vec2d &gl_l);
  explicit Domain_2(const Vec2d &gl_l, const Vec2i &gl_size,
                    const Vec2i &gl_cells_size, Vec2i &cells_size,
                    Vec2d &origin, Vec2b &flag_c);
  explicit Domain_2(const Vec2d &gl_l, const Vec2i &gl_cells_size,
                    Vec2i &cells_size, Vec2d &origin, Vec2b &flag_c);
  template <typename TNode>
  explicit Domain_2(const Vec2d &gl_l, CellListNode_2<TNode> **cl, double r_cut = 1.);

  void find_neighbor(Vec2i &rank, Vec2b &flag_comm, int neighbor[2][2]) const;

  void set_l(const Vec2i &gl_cells_size, Vec2i &cells_size,
             Vec2d &l, Vec2d &origin) const;

  void set_max_buf_size(int gl_par_num, double amplification);



  const Vec2d & l() const { return l_; }
  const Vec2d & gl_l() const { return gl_l_; }
  const Vec2d & gl_half_l() const { return gl_half_l_; }
  const Vec2d & origin() const { return origin_; }
  const Vec2b & flag_comm() const { return flag_comm_; }
  const Vec2i & domain_sizes() const { return gl_size_; }
  const Vec2i & domain_rank() const { return rank_; }
  const Vec2i & gl_cells_size() const { return gl_cells_size_; }
  const Vec2i & cells_size() const { return cells_size_; }
  const int & max_buf_size() const { return max_buf_size_; }

  static Vec_2<int> partition(const Vec2d &l, int n_proc);
  static Vec_2<int> partition(const Vec2d &l);

  int neighbor[2][2]{};
  Vec_2<block_t> inner_shell[2]{};
  Vec_2<block_t> outer_shell[2]{};
protected:
  Vec2d l_;
  Vec2d origin_{};
  Vec2d gl_l_;
  Vec2d gl_half_l_;

  Vec2i gl_size_{ 1, 1};
  Vec2i rank_{};
  Vec2b flag_comm_{};

  Vec2i gl_cells_size_{};
  Vec2i cells_size_{};

  int max_buf_size_ = 0;
};

template <typename TNode>
Domain_2::Domain_2(const Vec2d& gl_l, CellListNode_2<TNode>** cl, double r_cut)
  : gl_l_(gl_l), gl_half_l_(gl_l * 0.5), gl_size_(partition(gl_l_)) {
  Vec_2<double> cell_len{};
  CellListNode_2<TNode>::partition(gl_l_, r_cut, gl_cells_size_, cell_len);
  find_neighbor(rank_, flag_comm_, neighbor);
  set_l(gl_cells_size_, cells_size_, l_, origin_);
  set_comm_block(cells_size_, flag_comm_, inner_shell, outer_shell);
  *cl = new CellListNode_2<TNode>(cells_size_, cell_len, gl_l_, origin_, flag_comm_);
}
