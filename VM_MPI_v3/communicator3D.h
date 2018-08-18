#pragma once
#include "mpi.h"
#include "vect.h"
#include "cellList3D.h"
#include <algorithm>

struct block_t {
  Vec_3<int> beg;
  Vec_3<int> end;
};

void test_comm_velocity();

void find_shell(const Vec_3<int> &n, const Vec_3<int> &thickness,
                Vec_3<block_t> shell[2]);

void set_comm_block(const Vec_3<int> &cells_size, const Vec_3<bool> &flag_comm,
                    Vec_3<block_t> inner_shell[2], Vec_3<block_t> outer_shell[2]);

template <typename TNode>
void pack_ghost_par(double *buf, int &buf_size,
                    CellListNode_3<TNode>& cl,
                    const block_t& block) {
  int pos = 0;
  auto f_copy = [&pos, buf](TNode **head) {
    if (*head) {
      TNode *cur_node = *head;
      do {
        cur_node->copy(buf, pos);
        cur_node = cur_node->next;
      } while (cur_node);
    }
  };
  cl.for_each_cell(f_copy, block.beg, block.end);
  buf_size = pos;
}

template <typename TNode>
void unpack_ghost_par(const double *buf, int buf_size,
                      CellListNode_3<TNode>& cl,
                      std::vector<TNode> &p_arr,
                      int &n_ghost) {
  Vec_3<double> p1(Vec_3<double>(buf[0], buf[1], buf[2]));
  const Vec_3<double> offset = cl.get_offset(p1);
  const int n_new = buf_size / 6;
  const int n0 = p_arr.size();
  for (int i = 0; i < n_new; i++) {
    int ip = i + n0;
    p_arr.emplace_back(&buf[i * 6]);
    p_arr[ip].pos += offset;
    cl.add_node(p_arr[ip]);
    n_ghost++;
  }
}

template <typename TNode>
void pack_leaving_par(const std::vector<TNode> &p_arr,
                      std::vector<int> &vacant_pos,
                      CellListNode_3<TNode>& cl,
                      const block_t &block,
                      double *buf, int &buf_size) {
  const TNode* p0 = &p_arr[0];
  int buf_pos = 0;
  auto f_copy = [&buf_pos, &vacant_pos, p0, buf](TNode **head) {
    if (*head) {
      TNode *cur_node = *head;
      do {
        cur_node->copy(buf, buf_pos);
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
                        CellListNode_3<TNode>& cl,
                        std::vector<TNode> &p_arr,
                        std::vector<int> &vacant_pos) {  //! should be sorted in desending order
  const Vec_3<double> offset = cl.get_offset(Vec_3<double>(buf[0], buf[1], buf[2]));
  for (int buf_pos = 0; buf_pos < buf_size; buf_pos += 6) {
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

template <typename TPack, typename TUnpack, typename TFunc>
void comm_par(int prev_proc, int next_proc, int tag_bw, int tag_fw,
              const block_t &prev_block, const block_t &next_block,
              int max_buf_size, TPack pack, TUnpack unpack, TFunc do_sth) {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];
  for (int i = 0; i < 4; i++) {
    buf[i] = new double[max_buf_size];
    buf_size[i] = max_buf_size;
  }

  //! transfer data backward
  MPI_Irecv(buf[0], buf_size[0], MPI_DOUBLE, next_proc, tag_bw, MPI_COMM_WORLD, &req[0]);
  pack(buf[1], buf_size[1], prev_block);
  MPI_Isend(buf[1], buf_size[1], MPI_DOUBLE, prev_proc, tag_bw, MPI_COMM_WORLD, &req[1]);
  //! tarnsfer data forward
  MPI_Irecv(buf[2], buf_size[2], MPI_DOUBLE, prev_proc, tag_fw, MPI_COMM_WORLD, &req[2]);
  pack(buf[3], buf_size[3], next_block);
  MPI_Isend(buf[3], buf_size[3], MPI_DOUBLE, next_proc, tag_fw, MPI_COMM_WORLD, &req[3]);


  //! do something while waiting
  do_sth();

  MPI_Wait(&req[0], &stat[0]);
  MPI_Get_count(&stat[0], MPI_DOUBLE, &buf_size[0]);
  unpack(buf[0], buf_size[0]);

  MPI_Wait(&req[2], &stat[2]);
  MPI_Get_count(&stat[2], MPI_DOUBLE, &buf_size[2]);
  unpack(buf[2], buf_size[2]);

  MPI_Wait(&req[1], &stat[1]);
  MPI_Wait(&req[3], &stat[3]);

  for (int i = 0; i < 4; i++) {
    delete[] buf[i];
  }
}

template <typename TNode>
void comm_par_before_interact(const int neighbor_proc[3][2],
                              const Vec_3<block_t> inner_shell[2],
                              int max_buf_size,
                              std::vector<TNode> &p_arr,
                              CellListNode_3<TNode> &cl,
                              int& n_ghost) {
  n_ghost = 0;

  auto pack = [&cl](double *buf, int &buf_size, const block_t& block) {
    pack_ghost_par(buf, buf_size, cl, block);
  };

  auto unpack = [&n_ghost, &cl, &p_arr](double *buf, int buf_size) {
    unpack_ghost_par(buf, buf_size, cl, p_arr, n_ghost);
  };

  for (int direction = 0; direction < 3; direction++) {
    if (cl.flag_ext()[direction]) {
      const int prev_proc = neighbor_proc[direction][0];
      const int next_proc = neighbor_proc[direction][1];
      const auto & prev_block = inner_shell[0][direction];
      const auto & next_block = inner_shell[1][direction];
      comm_par(prev_proc, next_proc, 13, 31, prev_block, next_block,
               max_buf_size, pack, unpack, []() {});
    }
  }
}

template <typename TNode>
void clear_ghost_after_interact(CellListNode_3<TNode> &cl,
                                const Vec_3<block_t> outer_shell[2],
                                std::vector<TNode> &p_arr, int n_ghost) {
  for (int dim = 0; dim < 3; dim++) {
    if (cl.flag_ext()[dim]) {
      cl.clear(outer_shell[0][dim].beg, outer_shell[0][dim].end);
      cl.clear(outer_shell[1][dim].beg, outer_shell[1][dim].end);
    }
  }
  for (int i = 0; i < n_ghost; i++) {
    p_arr.pop_back();
  }
}

template <typename TPar>
void make_compact(std::vector<BiNode<TPar>> &p_arr,
                  std::vector<int> &vacancy,          //! in desending order
                  CellListNode_3<BiNode<TPar>> &cl) {
  int k = 0;
  while (k < vacancy.size()) {
    if (p_arr.size() - 1 == vacancy[k]) {
      p_arr.pop_back();
      k++;
    } else {
      auto i = vacancy.back();
      cl.replace(p_arr[i], p_arr.back());
      p_arr.pop_back();
      vacancy.pop_back();
    }
  }
}

template <typename T>
void make_compact(std::vector<T> &arr,
                  std::vector<int> &vacancy) {
  int k = 0;
  while (k < vacancy.size()) {
    if (arr.size() - 1 == vacancy[k]) {
      arr.pop_back();
      k++;
    } else {
      auto i = vacancy.back();
      arr[i] = arr.back();
      arr.pop_back();
      vacancy.pop_back();
    }
  }
}

template <typename TNode>
void comm_par_after_move(const int neighbor_proc[3][2],
                         const Vec_3<block_t> outer_shell[2],
                         int max_buf_size,
                         std::vector<TNode> &p_arr,
                         CellListNode_3<TNode>& cl) {
  std::vector<int> vacant_pos;
  vacant_pos.reserve(max_buf_size);
  auto pack = [&p_arr, &vacant_pos, &cl](double *buf, int &buf_size, const block_t& block) {
    pack_leaving_par(p_arr, vacant_pos, cl, block, buf, buf_size);
  };

  auto unpack = [&p_arr, &vacant_pos, &cl](double *buf, int buf_size) {
    unpack_arrived_par(buf, buf_size, cl, p_arr, vacant_pos);
  };

  auto sort_desending = [&vacant_pos]() {
    std::sort(vacant_pos.begin(), vacant_pos.end(), std::greater<int>());
  };

  for (int direction = 0; direction < 3; direction++) {
    if (cl.flag_ext()[direction]) {
      const int prev_proc = neighbor_proc[direction][0];
      const int next_proc = neighbor_proc[direction][1];
      const auto & prev_block = outer_shell[0][direction];
      const auto & next_block = outer_shell[1][direction];
      comm_par(prev_proc, next_proc, 24, 42, prev_block, next_block,
               max_buf_size, pack, unpack, sort_desending);
    }
  }

  make_compact(p_arr, vacant_pos, cl);
}