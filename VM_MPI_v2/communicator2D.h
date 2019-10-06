#pragma once
#include "mpi.h"
#include "vect.h"
#include "cellList2D.h"
#include <algorithm>

struct block_t {
  Vec_2<int> beg;
  Vec_2<int> end;
};

void find_shell(const Vec_2<int> &n, const Vec_2<int> &thickness, Vec_2<block_t> shell[2]);

Vec_2<int> decompose_domain(const Vec_2<double> &gl_l);

template <typename TNode>
void pack_ghost_par(double *buf, int &buf_size,
                    CellListNode_2<TNode>& cl,
                    const block_t& block) {
  int pos = 0;
  auto f_copy = [&pos, buf](TNode **head) {
    if (*head) {
      TNode *cur_node = *head;
      do {
        cur_node->copy_to(buf, pos);
        cur_node = cur_node->next;
      } while (cur_node);
    }
  };
  cl.for_each_cell(f_copy, block.beg, block.end);
  buf_size = pos;
}

template <typename TNode>
void unpack_ghost_par(const double *buf, int buf_size,
                      CellListNode_2<TNode>& cl,
                      std::vector<TNode> &p_arr,
                      int &n_ghost) {
  const Vec_2<double> offset = get_offset(Vec_2<double>(buf[0], buf[1]), cl);
  for (int buf_pos = 0; buf_pos < buf_size; buf_pos += 4) {
    auto idx_last = p_arr.size();
    p_arr.emplace_back(&buf[buf_pos]);
    p_arr[idx_last].pos += offset;
    cl.add_node(p_arr[idx_last]);
    n_ghost++;
  }
}

template <typename TNode>
void pack_leaving_par(const std::vector<TNode> &p_arr,
                      std::vector<int> &vacant_pos,
                      CellListNode_2<TNode>& cl,
                      const block_t &block,
                      double *buf, int &buf_size) {
  const TNode* p0 = &p_arr[0];
  int buf_pos = 0;
  auto f_copy = [&buf_pos, &vacant_pos, p0, buf](TNode **head) {
    if (*head) {
      TNode *cur_node = *head;
      do {
        cur_node->copy_to(buf, buf_pos);
        vacant_pos.push_back(cur_node - p0);
        cur_node = cur_node->next;
      } while (cur_node);
      *head = nullptr;
    }
  };
  cl.for_each_cell(f_copy, block.beg, block.end);
  buf_size = buf_pos;
}

template <typename TNode>
void unpack_arrived_par(const double *buf, int buf_size,
                        CellListNode_2<TNode>& cl,
                        std::vector<TNode> &p_arr,
                        std::vector<int> &vacant_pos) {  //! should be sorted in descending order
  const Vec_2<double> offset = get_offset(Vec_2<double>(buf[0], buf[1]), cl);
  for (int buf_pos = 0; buf_pos < buf_size; buf_pos += 4) {
    int idx;
    if (vacant_pos.empty()) {
      idx = p_arr.size();
      p_arr.emplace_back(&buf[buf_pos]);
    } else {
      idx = vacant_pos.back();
      p_arr[idx] = TNode(&buf[buf_pos]);
      vacant_pos.pop_back();
    }
    p_arr[idx].pos += offset;
    cl.add_node(p_arr[idx]);
  }
}

template <typename TNode, typename T1, typename T2>
void unpack_arrived_par(const double *buf, int buf_size,
                        CellListNode_2<TNode>& cl,
                        std::vector<TNode> &p_arr,
                        std::vector<int> &vacant_pos, //! should be sorted in descending order
                        std::vector<T1> &n_arr,
                        std::vector<Vec_2<T2>> &v_arr) {  
  const Vec_2<double> offset = get_offset(Vec_2<double>(buf[0], buf[1]), cl);
  for (int buf_pos = 0; buf_pos < buf_size; buf_pos += 4) {
    int idx;
    if (vacant_pos.empty()) {
      idx = p_arr.size();
      p_arr.emplace_back(&buf[buf_pos]);
    } else {
      idx = vacant_pos.back();
      p_arr[idx] = TNode(&buf[buf_pos]);
      vacant_pos.pop_back();
    }
    p_arr[idx].pos += offset;
    cl.add_node(p_arr[idx], n_arr, v_arr);
  }
}


class Communicator {
public:
  Communicator(const Vec_2<double> &gl_l, double rho0, double amplification,
               const Vec_2<int> &cells_size, const Vec_2<int> &domains_size);

  void find_neighbor(const Vec_2<int> &domain_rank,
                     const Vec_2<int> &domain_size,
                     Vec_2<bool> &flag_comm);

  int get_max_buf_size(double rho0, double amplification,
                        const Vec_2<double> &l) const;

  void set_comm_shell(const Vec_2<int> &cells_size);

  template <typename TPack, typename TUnpack, typename TFunc>
  void exchange_particle(int prev_proc, int next_proc, int tag_bw, int tag_fw,
                         const block_t &prev_block, const block_t &next_block,
                         TPack pack, TUnpack unpack, TFunc do_sth);

  template <typename TNode>
  void comm_before_cal_force(std::vector<TNode> &p_arr, CellListNode_2<TNode> &cl, int& n_ghost);

  template <typename TNode>
  void clear_padded_particles(CellListNode_2<TNode> &cl, std::vector<TNode> &p_arr, int n_ghost);

  template <typename TNode>
  void comm_after_integration(std::vector<TNode> &p_arr, CellListNode_2<TNode>& cl);

  template <typename TNode, typename T1, typename T2>
  void comm_after_integration(std::vector<TNode> &p_arr,
                              CellListNode_2<TNode>& cl,
                              std::vector<T1> &n_arr,
                              std::vector<Vec_2<T2>> &v_arr);

private:
  int tot_proc_ = 1;
  int my_rank_ = 0;

  Vec_2<bool> flag_comm_{};
  int neighbor_[2][2]{};
  Vec_2<block_t> inner_shell_[2]{};
  Vec_2<block_t> outer_shell_[2]{};

  double *buf_[4]{};
  int buf_size_[4]{};
  int max_buf_size_ = 0;

  std::vector<int> vacant_pos_;
};

template <typename TPack, typename TUnpack, typename TFunc>
void Communicator::exchange_particle(int prev_proc, int next_proc, int tag_bw, int tag_fw,
                                     const block_t& prev_block, const block_t& next_block,
                                     TPack pack, TUnpack unpack, TFunc do_sth) {
  MPI_Request req[4];
  MPI_Status stat[4];
  for (int i = 0; i < 4; i++) {
    buf_size_[i] = max_buf_size_;
  }

  //! transfer data backward
  MPI_Irecv(buf_[0], buf_size_[0], MPI_DOUBLE, next_proc, tag_bw, MPI_COMM_WORLD, &req[0]);
  pack(buf_[1], buf_size_[1], prev_block);
  MPI_Isend(buf_[1], buf_size_[1], MPI_DOUBLE, prev_proc, tag_bw, MPI_COMM_WORLD, &req[1]);
  //! transfer data forward
  MPI_Irecv(buf_[2], buf_size_[2], MPI_DOUBLE, prev_proc, tag_fw, MPI_COMM_WORLD, &req[2]);
  pack(buf_[3], buf_size_[3], next_block);
  MPI_Isend(buf_[3], buf_size_[3], MPI_DOUBLE, next_proc, tag_fw, MPI_COMM_WORLD, &req[3]);

  //! do something while waiting
  do_sth();

  MPI_Wait(&req[0], &stat[0]);
  MPI_Get_count(&stat[0], MPI_DOUBLE, &buf_size_[0]);
  unpack(buf_[0], buf_size_[0]);

  MPI_Wait(&req[2], &stat[2]);
  MPI_Get_count(&stat[2], MPI_DOUBLE, &buf_size_[2]);
  unpack(buf_[2], buf_size_[2]);

  MPI_Wait(&req[1], &stat[1]);
  MPI_Wait(&req[3], &stat[3]);
}

template <typename TNode>
void Communicator::comm_before_cal_force(std::vector<TNode>& p_arr,
                                         CellListNode_2<TNode>& cl, int& n_ghost) {
  n_ghost = 0;
  auto pack = [&cl](double *buf, int &buf_size, const block_t& block) {
    pack_ghost_par(buf, buf_size, cl, block);
  };

  auto unpack = [&n_ghost, &cl, &p_arr](double *buf, int buf_size) {
    int new_size = buf_size / 4 + p_arr.size();
    if (new_size > p_arr.capacity()) {
      cl.reserve_particles(p_arr, new_size);
    }
    unpack_ghost_par(buf, buf_size, cl, p_arr, n_ghost);
  };

  for (int direction = 0; direction < 2; direction++) {
    if (cl.flag_ext()[direction]) {
      const int prev_proc = neighbor_[direction][0];
      const int next_proc = neighbor_[direction][1];
      const auto & prev_block = inner_shell_[0][direction];
      const auto & next_block = inner_shell_[1][direction];
      exchange_particle(prev_proc, next_proc, 13, 31,
                        prev_block, next_block, pack, unpack, []() {});
    }
  }
}

template <typename TNode>
void Communicator::clear_padded_particles(CellListNode_2<TNode>& cl,
                                          std::vector<TNode>& p_arr, int n_ghost) {
  for (int dim = 0; dim < 2; dim++) {
    if (cl.flag_ext()[dim]) {
      cl.clear(outer_shell_[0][dim].beg, outer_shell_[0][dim].end);
      cl.clear(outer_shell_[1][dim].beg, outer_shell_[1][dim].end);
    }
  }
  for (int i = 0; i < n_ghost; i++) {
    p_arr.pop_back();
  }
}

template <typename TNode>
void Communicator::comm_after_integration(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl) {
  auto pack = [&p_arr, this, &cl](double *buf, int &buf_size, const block_t& block) {
    pack_leaving_par(p_arr, vacant_pos_, cl, block, buf, buf_size);
  };

  auto unpack = [&p_arr, this, &cl](double *buf, int buf_size) {
    int new_size = buf_size / 4 + p_arr.size() - vacant_pos_.size();
    if (new_size > p_arr.capacity()) {
      cl.reserve_particles(p_arr, new_size);
    }
    unpack_arrived_par(buf, buf_size, cl, p_arr, vacant_pos_);
  };

  auto sort_descending = [this]() {
    std::sort(vacant_pos_.begin(), vacant_pos_.end(), std::greater<int>());
  };

  for (int direction = 0; direction < 2; direction++) {
    if (cl.flag_ext()[direction]) {
      const int prev_proc = neighbor_[direction][0];
      const int next_proc = neighbor_[direction][1];
      const auto & prev_block = outer_shell_[0][direction];
      const auto & next_block = outer_shell_[1][direction];
      exchange_particle(prev_proc, next_proc, 24, 42, prev_block, next_block,
                        pack, unpack, sort_descending);
    }
  }

  cl.make_compact(p_arr, vacant_pos_);
}

template <typename TNode, typename T1, typename T2>
void Communicator::comm_after_integration(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl,
                                          std::vector<T1>& n_arr, std::vector<Vec_2<T2>>& v_arr) {
  auto pack = [&p_arr, this, &cl](double *buf, int &buf_size, const block_t& block) {
    pack_leaving_par(p_arr, vacant_pos_, cl, block, buf, buf_size);
  };

  auto unpack = [&p_arr, this, &cl, &n_arr, &v_arr](double *buf, int buf_size) {
    int new_size = buf_size / 4 + p_arr.size() - vacant_pos_.size();
    if (new_size > p_arr.capacity()) {
      cl.reserve_particles(p_arr, new_size);
    }
    unpack_arrived_par(buf, buf_size, cl, p_arr, vacant_pos_, n_arr, v_arr);
  };

  auto sort_descending = [this]() {
    std::sort(vacant_pos_.begin(), vacant_pos_.end(), std::greater<int>());
  };

  for (int direction = 0; direction < 2; direction++) {
    if (cl.flag_ext()[direction]) {
      const int prev_proc = neighbor_[direction][0];
      const int next_proc = neighbor_[direction][1];
      const auto & prev_block = outer_shell_[0][direction];
      const auto & next_block = outer_shell_[1][direction];
      exchange_particle(prev_proc, next_proc, 24, 42, prev_block, next_block,
        pack, unpack, sort_descending);
    }
  }
  cl.make_compact(p_arr, vacant_pos_);
}
