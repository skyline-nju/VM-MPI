#pragma once
#include "vect.h"
#include "comn.h"

class Domain_2 {
public:
  explicit Domain_2(const Vec_2<double> &gl_l, const Vec_2<int> &domains_size);

  void tangle(Vec_2<double> &pos) const;

  void untangle(Vec_2<double>& v) const;

  const Vec_2<double> & l() const { return l_; }
  const Vec_2<double> & gl_l() const { return gl_l_; }
  const Vec_2<double> & gl_half_l() const { return gl_half_l_; }
  const Vec_2<double> & origin() const { return origin_; }
  const Vec_2<bool> & flag_comm() const { return flag_comm_; }
  const Vec_2<int> &rank() const { return rank_; }
protected:
  Vec_2<double> l_{};
  Vec_2<double> origin_{};
  Vec_2<double> gl_l_;
  Vec_2<double> gl_half_l_;
  Vec_2<bool> flag_comm_{};
  Vec_2<int> rank_{};
};


inline void Domain_2::tangle(Vec_2<double>& pos) const {
  if (!flag_comm_.x) {
    tangle_1(pos.x, gl_l_.x);
  }
  if (!flag_comm_.y) {
    tangle_1(pos.y, gl_l_.y);
  }
}

inline void Domain_2::untangle(Vec_2<double>& v) const {
  if (!flag_comm_.x) {
    untangle_1(v.x, gl_l_.x, gl_half_l_.x);
  }
  if (flag_comm_.y) {
    untangle_1(v.y, gl_l_.y, gl_half_l_.y);
  }
}
