#pragma once
#include "mpi.h"
#include "vect.h"
#include "cellList3D.h"
#include <algorithm>

struct block_t {
  Vec_3<int> beg;
  Vec_3<int> end;
};

int get_proc_num();

int get_proc_rank();

void find_neighbor_proc(const Vec_3<int> &domain_size,
                        Vec_3<int> &domain_rank,
                        int neighbor_proc[3][2],
                        Vec_3<bool> &flag_ext,
                        int test_rank = -1);

void divide_domain(const Vec_3<int>& domain_size,
                   const Vec_3<int>& domain_rank,
                   const Vec_3<int>& gl_cell_size,
                   Vec_3<int>& my_cell_size,
                   Vec_3<int>& my_first_cell);

int divide_par_evenly(int gl_np);

void set_comm_block(const Vec_3<int> &cell_size, const Vec_3<bool> &flag_comm);

const Vec_3<block_t>* get_inner_block();

const Vec_3<block_t>* get_outer_block();

template <typename TNode>
void pack_ghost_par(double *buf, int &buf_size,
                    const CellListNode_3<TNode>& cl,
                    const block_t& block) {
  int pos = 0;
  auto f_copy = [&pos](const TNode *head) {
    TNode *cur_node = head;
    while(cur_node) {
      buf[pos]     = cur_node->pos.x;
      buf[pos + 1] = cur_node->pos.y;
      buf[pos + 2] = cur_node->pos.z;
      buf[pos + 3] = cur_node->ori.x;
      buf[pos + 4] = cur_node->ori.y;
      buf[pos + 5] = cur_node->ori.z;
      pos += 6;
      cur_node = cur_node->next;
    }
  };
  cl.for_each_node(f_copy, block.beg, block.end);
  buf_size = pos;
}

template <typename TNode>
void unpack_ghost_par(const double *buf, int buf_size,
                      CellListNode_3<TNode>& cl,
                      std::vector<TNode> &p_arr,
                      int &n_ghost) {
  const Vec_3<double> offset = cl.get_offset(Vec_3<double>(buf[0], buf[1], buf[2]));
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
void pack_leave_par(const std::vector<TNode> &p_arr,
                    std::vector<int> &ip_out,
                    const CellListNode_3<TNode>& cl,
                    const block_t &block,
                    double *buf, int &buf_size) {
  TNode* p0 = &p_arr[0];
  int pos = 0;
  auto f_copy = [&pos](const TNode *head) {
    TNode *cur_node = head;
    while (cur_node) {
      buf[pos]     = cur_node->pos.x;
      buf[pos + 1] = cur_node->pos.y;
      buf[pos + 2] = cur_node->pos.z;
      buf[pos + 3] = cur_node->ori.x;
      buf[pos + 4] = cur_node->ori.y;
      buf[pos + 5] = cur_node->ori.z;
      pos += 6;
      ip_out.push_back(cur_node - p0);
      cur_node = cur_node->next;
    }
  };
  cl.for_each_node_clean(f_copy, block.beg, block.end);
  buf_size = pos;
}

template <typename TNode>
void unpack_arrive_par(const double *buf, int buf_size,
                       CellListNode_3<TNode>& cl,
                       std::vector<TNode> &p_arr,
                       std::vector<int> &ip_out_decending) {
  const Vec_3<double> offset = cl.get_offset(Vec_3<double>(buf[0], buf[1], buf[2]));
  const int n_new = buf_size / 6;
  for (int i = 0; i < n_new; i++) {
    int idx_new;
    if (ip_out_decending.empty()) {
      idx_new = p_arr.size();
      p_arr.emplace_back(&buf[i * 6]);
    } else {
      idx_new = ip_out_decending.back();
      p_arr[idx_new] = TNode(&buf[i * 6]);
      ip_out_decending.pop_back();
    }
    p_arr[idx_new].pos += offset;
    cl.add_node(p_arr[idx_new]);
  }
}

template <typename TPack, typename TUnpack, typename TFunc>
void par_comm(int direction, int neighbor_proc[3][2],
              const Vec_3<block_t> block[2], int max_buf_size,
              TPack pack, TUnpack unpack, TFunc do_something) {
  int my_rank = get_proc_rank();
  const int prev_rank = neighbor_proc[direction][0];
  const int next_rank = neighbor_proc[direction][1];
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];
  for (int i = 0; i < 4; i++) {
    buf[i] = new double[max_buf_size];
    buf_size[i] = max_buf_size;
  }

  MPI_Irecv(buf[0], buf_size[0], MPI_DOUBLE, next_rank,
            32, MPI_COMM_WORLD, &req[0]);
  pack(buf[1], &buf_size[1], block[0][direction]);
  MPI_Isend(buf[1], buf_size[1], MPI_DOUBLE, prev_rank,
            32, MPI_COMM_WORLD, &req[1]);

  MPI_Irecv(buf[2], buf_size[2], MPI_DOUBLE, prev_rank,
            23, MPI_COMM_WORLD, &req[2]);
  pack(buf[3], &buf_size[3], block[1][direction]);
  MPI_Isend(buf[3], buf_size[3], MPI_DOUBLE, next_rank,
            23, MPI_COMM_WORLD, &req[3]);

  do_something();

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
void par_comm_before_interact(int neighbor_proc[3][2],
                              int max_buf_size,
                              std::vector<TNode> &p_arr,
                              CellListNode_3<TNode> &cl,
                              int& n_ghost) {
  auto pack = [&cl](double *buf, int &buf_size, const block_t& block) {
    pack_ghost_par(buf, buf_size, cl, block);
  };

  auto unpack = [&n_ghost, &cl, &p_arr](double *buf, int buf_size) {
    unpack_ghost_par(buf, buf_size, cl, p_arr, n_ghost);
  };

  n_ghost = 0;
  const Vec_3<block_t> *inner_block = get_inner_block();
  for (int direction = 0; direction < 3; direction++) {
    if (cl.flag_ext()[direction])
      par_comm(direction, neighbor_proc, inner_block, max_buf_size,
               pack, unpack, [](){});
  }
}

template <typename TNode>
void clear_after_interact(CellListNode_3<TNode> &cl,
  std::vector<TNode> &p_arr, int n_ghost) {
  const Vec_3<block_t> *outer_block = get_outer_block();
  for (int dim = 0; dim < 3; dim++) {
    if (cl.flag_ext()[dim]) {
      cl.clear(outer_block[0][dim].beg, outer_block[0][dim].end);
      cl.clear(outer_block[1][dim].beg, outer_block[1][dim].end);
    }
  }
  for (int i = 0; i < n_ghost; i++) {
    p_arr.pop_back();
  }
}


template <typename TPar>
void make_compact(std::vector<BiNode<TPar>> &p_arr,
                  std::vector<int> &vacant_pos,          //! in desending order
                  CellListNode_3<BiNode<TPar>> &cl) {
  int k = 0;
  while(k < vacant_pos.size()) {
    while (p_arr.size() - 1 == vacant_pos[k]) {
      p_arr.pop_back();
      k++;
    }
    auto first_vac_pos = vacant_pos.back();
    cl.replace(p_arr[first_vac_pos], p_arr.back());
    p_arr.pop_back();
    vacant_pos.pop_back();
  }
}

template <typename T>
void make_compact(std::vector<T> &arr,
                  std::vector<int> &vacant_pos) {
  int k = 0;
  while (k < vacant_pos.size()) {
    while (arr.size() - 1 == vacant_pos[k]) {
      arr.pop_back();
      k++;
    }
    auto first_vac_pos = vacant_pos.back();
    arr[first_vac_pos] = arr.back();
    arr.pop_back();
    vacant_pos.pop_back();
  }
}

template <typename TNode>
void par_comm_after_move(int neighbor_proc[3][2],
                         int max_buf_size,
                         std::vector<TNode> &p_arr,
                         CellListNode_3<TNode>& cl) {
  std::vector<int> vacant_pos;
  vacant_pos.reserve(max_buf_size * 6);
  auto pack = [&p_arr, &vacant_pos, &cl](double *buf, int &buf_size, const block_t& block) {
    pack_leave_par(p_arr, vacant_pos, cl, block, buf, buf_size);
  };

  auto unpack = [&p_arr, &vacant_pos, &cl](double *buf, int buf_size) {
    unpack_arrive_par(buf, buf_size, cl, p_arr, vacant_pos);
  };

  auto sort_desending = [&vacant_pos]() {
    std::sort(vacant_pos.begin(), vacant_pos.end(), std::greater<int>());
  };
  const Vec_3<block_t> *outer_block = get_outer_block();
  for (int direction = 2; direction > -1; direction--) {
    if (cl.flag_ext()[direction])
      comm(direction, neighbor_proc, outer_block, max_buf_size,
           pack, unpack, sort_desending);
  }
  make_compact(p_arr, vacant_pos, cl);
}