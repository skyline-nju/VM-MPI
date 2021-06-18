#pragma once
#ifdef USE_MPI

#include "mpi.h"
#include "vect.h"
#include "cellList2D.h"
#include <algorithm>

template <typename T>  // default type for T is int
void find_shell(const Vec_2<T>& n, const Vec_2<T>& thickness, Vec_2<RectBlock_2<T>> shell[2]) {
  for (int ori = 0; ori < 2; ori++) {
    if (thickness[ori]) {
      for (int dim = 0; dim < 2; dim++) {
        if (dim == ori) {
          shell[0][ori].beg[dim] = 0;
          shell[0][ori].end[dim] = 1;
          shell[1][ori].beg[dim] = n[dim] - 1;
          shell[1][ori].end[dim] = n[dim];
        } else if (dim > ori) {
          shell[0][ori].beg[dim] = shell[1][ori].beg[dim] = 0;
          shell[0][ori].end[dim] = shell[1][ori].end[dim] = n[dim];
        } else {
          shell[0][ori].beg[dim] = shell[1][ori].beg[dim] = thickness[dim];
          shell[0][ori].end[dim] = shell[1][ori].end[dim] = n[dim] - thickness[dim];
        }
      }
    }
  }
}

template <typename TNode>
void pack_ghost_par(double *buf, int &buf_size,
                    CellListNode_2<TNode>& cl,
                    const RectBlock_2<int>& block) {
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
  //const Vec_2<double> offset = get_offset(Vec_2<double>(buf[0], buf[1]), cl);
  const Vec_2<double> offset = cl.get_pos_offset(Vec_2<double>(buf[0], buf[1]));
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
                      const RectBlock_2<int> &block,
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
  const Vec_2<double> offset = cl.get_pos_offset(Vec_2<double>(buf[0], buf[1]));
  
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
  const Vec_2<double> offset = cl.get_pos_offset(Vec_2<double>(buf[0], buf[1]));

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


class Communicator_2 {
public:
  template <class TDomain>
  Communicator_2(const TDomain& dm, double rho0, double amplification);

  template <class TDomain, class TGrid>
  Communicator_2(const TDomain& dm, const TGrid& grid, double rho0, double amplification);

  template <typename T>
  int get_max_buf_size(const T rho0, double amplification,
                       const Vec_2<double> &l) const;

  template <typename T>
  void set_comm_shell(const Vec_2<T> &cells_size);

  template <typename TPack, typename TUnpack, typename TFunc>
  void exchange_particle(int prev_proc, int next_proc, int tag_bw, int tag_fw,
                         const RectBlock_2<int> &prev_block, const RectBlock_2<int> &next_block,
                         TPack pack, TUnpack unpack, TFunc do_sth);

  template <typename TNode>
  void comm_before_cal_force(std::vector<TNode> &p_arr, CellListNode_2<TNode> &cl, int& n_ghost);

  template <typename TNode>
  void clear_padded_particles(CellListNode_2<TNode> &cl, std::vector<TNode> &p_arr, int n_ghost);

  template <typename TNode>
  void comm_after_integration(std::vector<TNode> &p_arr, CellListNode_2<TNode>& cl);

private:
  int tot_proc_ = 1;
  int my_rank_ = 0;

  Vec_2<bool> flag_comm_{};
  int neighbor_[2][2]{};
  Vec_2<RectBlock_2<int>> inner_shell_[2]{};
  Vec_2<RectBlock_2<int>> outer_shell_[2]{};

  double *buf_[4]{};
  int buf_size_[4]{};
  int max_buf_size_ = 0;

  std::vector<int> vacant_pos_;
  MPI_Comm comm_;
};

template <typename T>
int Communicator_2::get_max_buf_size(const T rho0, double amplification,
                                     const Vec_2<double>& l) const {
  std::vector<double> area;
  if (flag_comm_.x) {
    area.push_back(l.y);
  }
  if (flag_comm_.y) {
    area.push_back(l.x);
  }

  int max_buf_size;
  if (area.empty()) {
    max_buf_size = 0;

  } else {
    std::sort(area.begin(), area.end(), [](double x, double y) {return x > y; });

    int n0 = int(rho0 * area[0] * amplification);
    max_buf_size = 4 * n0;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
      std::cout << "max area = " << area[0] << std::endl;
      std::cout << "max particle number per communication: " << n0 << " particles" << std::endl;
    }
  }
  return max_buf_size;
}

template <typename T>
void Communicator_2::set_comm_shell(const Vec_2<T>& cells_size) {
  Vec_2<int> thickness{};
  for (int dim = 0; dim < 2; dim++) {
    thickness[dim] = flag_comm_[dim] ? 1 : 0;
  }

  find_shell(cells_size, -thickness, inner_shell_);
  for (int ori = 0; ori < 2; ori++) {
    inner_shell_[0][ori].beg += thickness;
    inner_shell_[0][ori].end += thickness;
    inner_shell_[1][ori].beg += thickness;
    inner_shell_[1][ori].end += thickness;
  }

  const Vec_2<int> extended_cells_size = cells_size + thickness * 2;
  find_shell(extended_cells_size, thickness, outer_shell_);
}

template<class TDomain>
Communicator_2::Communicator_2(const TDomain& dm, double rho0, double amplification):
                               flag_comm_(dm.flag_comm()), comm_(dm.comm()) {
  my_rank_ = dm.rank().x + dm.rank().y * dm.size().x;
  tot_proc_ = dm.size().x * dm.size().y;

  dm.find_neighbor(neighbor_);
  set_comm_shell(dm.grid().n());
  max_buf_size_ = get_max_buf_size(rho0, amplification, dm.l());
  for (int i = 0; i < 4; i++) {
    buf_[i] = new double[max_buf_size_];
    buf_size_[i] = max_buf_size_;
  }

  vacant_pos_.reserve(max_buf_size_);
}

template <class TDomain, class TGrid>
Communicator_2::Communicator_2(const TDomain& dm, const TGrid& grid, double rho0, double amplification):
                              flag_comm_(dm.proc_size().x > 1, dm.proc_size().y > 1), comm_(dm.comm()) {
  my_rank_ = dm.proc_rank().x + dm.proc_rank().y * dm.proc_size().x;
  tot_proc_ = dm.proc_size().x * dm.proc_size().y;

  dm.find_neighbor(neighbor_);
  set_comm_shell(grid.n());
  max_buf_size_ = get_max_buf_size(rho0, amplification, dm.l());
  for (int i = 0; i < 4; i++) {
    buf_[i] = new double[max_buf_size_];
    buf_size_[i] = max_buf_size_;
  }

  vacant_pos_.reserve(max_buf_size_);
}
template <typename TPack, typename TUnpack, typename TFunc>
void Communicator_2::exchange_particle(int prev_proc, int next_proc, int tag_bw, int tag_fw,
                                       const RectBlock_2<int>& prev_block, const RectBlock_2<int>& next_block,
                                       TPack pack, TUnpack unpack, TFunc do_sth) {
  MPI_Request req[4];
  MPI_Status stat[4];
  for (int i = 0; i < 4; i++) {
    buf_size_[i] = max_buf_size_;
  }

  //! transfer data backward
  MPI_Irecv(buf_[0], buf_size_[0], MPI_DOUBLE, next_proc, tag_bw, comm_, &req[0]);
  pack(buf_[1], buf_size_[1], prev_block);
  MPI_Isend(buf_[1], buf_size_[1], MPI_DOUBLE, prev_proc, tag_bw, comm_, &req[1]);
  //! transfer data forward
  MPI_Irecv(buf_[2], buf_size_[2], MPI_DOUBLE, prev_proc, tag_fw, comm_, &req[2]);
  pack(buf_[3], buf_size_[3], next_block);
  MPI_Isend(buf_[3], buf_size_[3], MPI_DOUBLE, next_proc, tag_fw, comm_, &req[3]);

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
void Communicator_2::comm_before_cal_force(std::vector<TNode>& p_arr,
                                           CellListNode_2<TNode>& cl, int& n_ghost) {
  n_ghost = 0;
  auto pack = [&cl](double *buf, int &buf_size, const RectBlock_2<int>& block) {
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
void Communicator_2::clear_padded_particles(CellListNode_2<TNode>& cl,
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
void Communicator_2::comm_after_integration(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl) {
  auto pack = [&p_arr, this, &cl](double *buf, int &buf_size, const RectBlock_2<int>& block) {
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

#endif