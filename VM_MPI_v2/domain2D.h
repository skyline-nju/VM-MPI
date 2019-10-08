#pragma once
#include "vect.h"
#include "comn.h"

Vec_2<int> get_proc_rank(const Vec_2<int>& proc_size);

class Grid_2 {
public:
  Grid_2(const Vec_2<double>& gl_l, const Vec_2<int>& proc_size, const Vec_2<int>& proc_rank, double r_cut);

  const Vec_2<int>& n() const { return n_; }
  const Vec_2<int>& gl_n() const { return gl_n_; }
  const Vec_2<int>& origin() const { return origin_; }
  const Vec_2<double>& lc() const { return lc_; }
private:
  Vec_2<int> n_;
  Vec_2<int> gl_n_;
  Vec_2<int> origin_;
  Vec_2<double> lc_; // length of each cell
};

class SubDomain_2 {
public:
  SubDomain_2(const Vec_2<double>& gl_l, const Vec_2<int>& proc_size, double r_cut);

  void tangle(Vec_2<double>& pos) const;
  void untangle(Vec_2<double>& v) const;

  void find_neighbor(int(*neighbor)[2]) const;

  const Vec_2<int>& size() const { return size_; }
  const Vec_2<int>& rank() const { return rank_; }
  const Grid_2& grid() const { return grid_; }
  const Vec_2<bool>& flag_comm() const { return flag_comm_; }

  const Vec_2<double>& gl_l() const { return gl_l_; }
  const Vec_2<double>& l() const { return l_; }
  const Vec_2<double>& origin() const { return origin_; }

  static const Vec_2<int> decompose(const Vec_2<double>& gl_l);
private:
  Vec_2<int> size_;
  Vec_2<int> rank_;
  Grid_2 grid_;
  Vec_2<bool> flag_comm_;
  Vec_2<double> gl_l_;
  Vec_2<double> half_gl_l_;
  Vec_2<double> l_;
  Vec_2<double> origin_;
};

inline void SubDomain_2::tangle(Vec_2<double>& pos) const {
  if (!flag_comm_.x) {
    tangle_1(pos.x, gl_l_.x);
  }
  if (!flag_comm_.y) {
    tangle_1(pos.y, gl_l_.y);
  }
}

inline void SubDomain_2::untangle(Vec_2<double>& v) const {
  if (!flag_comm_.x) {
    untangle_1(v.x, gl_l_.x, half_gl_l_.x);
  }
  if (flag_comm_.y) {
    untangle_1(v.y, gl_l_.y, half_gl_l_.y);
  }
}

