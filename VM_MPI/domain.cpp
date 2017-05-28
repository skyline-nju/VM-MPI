#include "domain.h"
#include <iostream>
#include <cmath>
using namespace std;

SubDomain::SubDomain(double Lx_domain, double Ly_domain, int nrows_subdomain,
                     int rank, unsigned long long seed){
  Lx = Lx_domain;
  Ly = Ly_domain;
  double dLy = Ly / nrows_subdomain;
  yl = rank * dLy - 1;
  yh = (rank + 1) * dLy + 1;
  ncols = int(Lx);
  nrows = int(dLy) + 2;
  int tot_cells = ncols * nrows;
  cell = new Cell[tot_cells];
  for (int i = 0; i < tot_cells; i++) {
    cell[i].find_neighbor(cell, i, ncols, nrows);
  }
  myran = new Ran(seed + rank);
  tot_rank = nrows_subdomain;
  myrank = rank;
  pre_rank = rank == 0 ? nrows_subdomain - 1 : rank - 1;
  next_rank = rank == nrows_subdomain - 1 ? 0 : rank + 1;
  cout << "subdomain " << rank << "\tpre rank " << pre_rank;
  cout << "\tnext rank " << next_rank << endl;
  cout << "\tncols = " << ncols << "\tnrows = " << nrows << endl;
}

void SubDomain::create_particle_random(int nPar) {
  double Ly0 = yl + 1;
  double Ly1 = yh - 1;
  particle.reserve(int(nPar * 3));
  cout << "max particles of rank " << myrank << " is " << particle.capacity() << endl;
  for (int i = 0; i < nPar; i++) {
    particle.push_back(Node());
    particle[i].x = myran->doub() * Lx;
    particle[i].y = myran->doub() * (Ly1 - Ly0) + Ly0;
    double theta = myran->doub() * 2 * PI;
    particle[i].vx = particle[i].vx0 = cos(theta);
    particle[i].vy = particle[i].vy0 = sin(theta);
    particle[i].par_idx = i;
  }
  create_cell_list();
  MAX_BUFF_SIZE = nPar / (nrows - 2) * 10 * 4;
}

void SubDomain::update_velocity_by_row(int row) {
  if (row > 0) {
    int j = row * ncols;
    { /* leftmost column */
      int i = j;
      if (cell[i].head) {
        cell[i].interact();
        cell[i].interact(cell[i].neighbor[0]);
        cell[i].interact(cell[i].neighbor[1], -Lx, 0);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3]);
      }
    }
    for (int col = 0; col < ncols - 1; col++) {
      int i = col + j;
      if (cell[i].head) {
        cell[i].interact();
        cell[i].interact(cell[i].neighbor[0]);
        cell[i].interact(cell[i].neighbor[1]);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3]);
      }
    }
    { /* rightmost row */
      int i = j + ncols - 1;
      if (cell[i].head) {
        cell[i].interact();
        cell[i].interact(cell[i].neighbor[0], Lx, 0);
        cell[i].interact(cell[i].neighbor[1]);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3], Lx, 0);
      }
    }
  } else {
    { /* col = 0 */
      int i = 0;
      if (cell[i].head) {
        cell[i].interact(cell[i].neighbor[1], -Lx, 0);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3]);
      }
    }
    for (int col = 1; col < ncols - 1; col++) {
      int i = col;
      if (cell[i].head) {
        cell[i].interact(cell[i].neighbor[1]);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3]);
      }
    }
    { /* col = ncols - 1 */
      int i = ncols - 1;
      if (cell[i].head) {
        cell[i].interact(cell[i].neighbor[1]);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3], Lx, 0);
      }
    }
  }
}

void SubDomain::update_velocity_inner_rows() {
  for (int row = 1; row < nrows - 2; row++) {
    update_velocity_by_row(row);
  }
}

void SubDomain::update_position_inner_rows(double eta) {
  for (int row = 1; row < nrows - 1; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      cell[col + j].head = nullptr;
      cell[col + j].size = 0;
    }
  }
  for (size_t i = 0, size = particle.size(); i < size; i++) {
    if (!particle[i].is_empty) {
      if (!particle[i].is_moved) {
        double noise = eta * 2 * PI * (myran->doub() - 0.5);
        particle[i].update_coor(noise, Lx, yl);
      }
      particle[i].is_moved = false;
      if (!particle[i].is_ghost)
        cell[particle[i].cell_idx].push_front(&particle[i]);
    }
  }
}

void SubDomain::update_position_edge_row(double eta, int row) {
  if (row != 1 && row != nrows - 2) {
    cout << "Error, row = " << row << " should be equal to "
         << 1 << " or " << nrows - 2 << endl;
    exit(1);
  }
  std::vector<Node *> ghost;
  ghost.reserve(MAX_BUFF_SIZE / 8);
  for (int i = row * ncols, n = ncols + row * ncols; i < n; i++) {
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        double noise = eta * 2 * PI * (myran->doub() - 0.5);
        curNode->update_coor(noise, Lx, yl);
        if (curNode->y < yl + 1 || curNode->y >= yh - 1) {
          curNode->is_ghost = true;
          ghost.push_back(curNode);
        }
        curNode = curNode->next;
      } while (curNode);
    }
  }
  for (int i = 0, size = ghost.size(); i < size; i++) {
    int c_idx = ghost[i]->cell_idx;
    cell[c_idx].push_front(ghost[i]);
  }
}

void SubDomain::create_cell_list() {
  for (size_t i = 0, size = particle.size(); i < size; i++) {
    particle[i].is_moved = false;
    particle[i].is_empty = false;
    particle[i].is_ghost = false;
    particle[i].new_arrival = false;
    particle[i].cal_cell_idx(yl, ncols);
    cell[particle[i].cell_idx].push_front(&particle[i]);
  }
}

void SubDomain::remove_ghost_particle(int row) {
  for (int col = 0; col < ncols; col++) {
    int idx = col + row * ncols;
    if (cell[idx].head) {
      Node *curNode = cell[idx].head;
      do {
        curNode->is_empty = true;
        empty_pos.push(curNode->par_idx);
        curNode = curNode->next;
      } while (curNode);
      cell[idx].head = nullptr;
      cell[idx].size = 0;
    }
  }
}


int SubDomain::count_valid_particle() {
  int count = 0;
  for (auto iter = particle.begin(); iter != particle.end(); ++iter) {
    if (!(*iter).is_empty && !(*iter).is_ghost)
      count += 1;
  }
  return count;
}

void SubDomain::pack(int row, double *buff, int &buff_size) {
  int pos = 0;
  double dy = 0;
  if (myrank == 0 && row <= 1) {
    dy = Ly;
  } else if (myrank == tot_rank - 1 && row >= nrows - 2) {
    dy = -Ly;
  }
  for (int col = 0, j = row * ncols; col < ncols; col++) {
    int i = col + j;
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        if (!curNode->new_arrival) {
          buff[pos++] = curNode->x;
          buff[pos++] = curNode->y + dy;
          buff[pos++] = curNode->vx;
          buff[pos++] = curNode->vy;
          if (curNode->y < -1 || curNode->y >= Ly + 1) {
            cout << "rank = " << myrank << "\ty = " << curNode->y << endl;
          }
        } else {
          curNode->new_arrival = false;
        }
        curNode = curNode->next;
        if (pos >= MAX_BUFF_SIZE) {
          cout << "Error, pos = " << pos << " is larger than MAX_BUFF_SIZE = " << MAX_BUFF_SIZE << endl;
          exit(1);
        }
      } while (curNode);
    }
  }
  buff_size = pos;
}

void SubDomain::unpack(int row, const double *buff, int buff_size) {
  int j = row * ncols;
  bool ghost;
  if (row == 0 || row == nrows - 1) {
    ghost = true;
  } else if (row == 1 || row == nrows -2) {
    ghost = false;
  } else {
    cout << "Error, wrong row to uppack\n";
    exit(1);
  }
  int nPar = buff_size / 4;
  for (unsigned int i = 0; i < nPar; i++) {
    double x = buff[i * 4];
    double y = buff[i * 4 + 1];
    double vx = buff[i * 4 + 2];
    double vy = buff[i * 4 + 3];
    int cell_idx = int(x) + j;
    if (int(y - yl) != row) {
      cout << "Error, rank = " << myrank << "\ty = " << y << "\tint(y - yl) = "
        << int(y - yl) << "\trow = " << row << endl;
      exit(1);
    }
    int par_idx;
    if (!empty_pos.empty()) {
      par_idx = empty_pos.top();
      empty_pos.pop();
    } else {
      par_idx = particle.size();
      particle.push_back(Node());
    }
    particle[par_idx].x = x;
    particle[par_idx].y = y;
    particle[par_idx].vx = particle[par_idx].vx0 = vx;
    particle[par_idx].vy = particle[par_idx].vy0 = vy;
    particle[par_idx].cell_idx = cell_idx;
    particle[par_idx].par_idx = par_idx;
    particle[par_idx].is_moved = false;
    particle[par_idx].is_empty = false;
    particle[par_idx].is_ghost = ghost;
    particle[par_idx].new_arrival = true;
    cell[cell_idx].push_front(&particle[par_idx]);
  }
}

void SubDomain::update_velocity_MPI() {
  MPI_Request reqs[4];
  MPI_Status stats[4];
  int tag_forward = 99;
  int tag_backward = 1;

  /* transfer data forward */
  double *buff_to_next = new double[MAX_BUFF_SIZE];
  double *buff_from_pre = new double[MAX_BUFF_SIZE];
  int size_to_next;
  pack(nrows - 2, buff_to_next, size_to_next);
  MPI_Irecv(buff_from_pre, MAX_BUFF_SIZE, MPI_DOUBLE, pre_rank,
            tag_forward, MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(buff_to_next, size_to_next, MPI_DOUBLE, next_rank,
            tag_forward, MPI_COMM_WORLD, &reqs[1]);

  /* transfer data backward */
  double *buff_to_pre = new double[MAX_BUFF_SIZE];
  double *buff_from_next = new double[MAX_BUFF_SIZE];
  int size_to_pre;
  pack(1, buff_to_pre, size_to_pre);
  MPI_Irecv(buff_from_next, MAX_BUFF_SIZE, MPI_DOUBLE, next_rank,
            tag_backward, MPI_COMM_WORLD, &reqs[2]);
  MPI_Isend(buff_to_pre, size_to_pre, MPI_DOUBLE, pre_rank,
            tag_backward, MPI_COMM_WORLD, &reqs[3]);

  /* update velocity of particles located in inner rows */
  update_velocity_inner_rows();

  /* unpack data to row = 0 and update velocity */
  MPI_Wait(&reqs[0], &stats[0]);
  int size_from_pre;
  MPI_Get_count(&stats[0], MPI_DOUBLE, &size_from_pre);
  unpack(0, buff_from_pre, size_from_pre);
  update_velocity_by_row(0);
  remove_ghost_particle(0);

  /* unpack data to row = nrows-1 and update velocity */
  MPI_Wait(&reqs[2], &stats[2]);
  int size_from_next;
  MPI_Get_count(&stats[2], MPI_DOUBLE, &size_from_next);
  unpack(nrows - 1, buff_from_next, size_from_next);
  update_velocity_by_row(nrows - 2);
  remove_ghost_particle(nrows - 1);

  /* wait until finishing sending */
  MPI_Wait(&reqs[1], &stats[1]);
  MPI_Wait(&reqs[3], &stats[3]);

  /* free memory */
  delete[] buff_from_next;
  delete[] buff_from_pre;
  delete[] buff_to_next;
  delete[] buff_to_pre;
}

void SubDomain::update_position_MPI(double eta) {
  MPI_Request reqs[4];
  MPI_Status stats[4];
  int tag_forward = 88;
  int tag_backward = 2;

  /* update position at row = nrows-2 and transfer data forward */
  update_position_edge_row(eta, nrows - 2);
  double *buff_to_next = new double[MAX_BUFF_SIZE];
  double *buff_from_pre = new double[MAX_BUFF_SIZE];
  int size_to_next;
  pack(nrows - 1, buff_to_next, size_to_next);
  MPI_Irecv(buff_from_pre, MAX_BUFF_SIZE, MPI_DOUBLE, pre_rank,
            tag_forward, MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(buff_to_next, size_to_next, MPI_DOUBLE, next_rank,
            tag_forward, MPI_COMM_WORLD, &reqs[1]);

  /* update positon at row = 1 and transfer data backward */
  update_position_edge_row(eta, 1);
  double *buff_to_pre = new double[MAX_BUFF_SIZE];
  double *buff_from_next = new double[MAX_BUFF_SIZE];
  int size_to_pre;
  pack(0, buff_to_pre, size_to_pre);
  MPI_Irecv(buff_from_next, MAX_BUFF_SIZE, MPI_DOUBLE, next_rank,
            tag_backward, MPI_COMM_WORLD, &reqs[2]);
  MPI_Isend(buff_to_pre, size_to_pre, MPI_DOUBLE, pre_rank,
            tag_backward, MPI_COMM_WORLD, &reqs[3]);

  /* update position of particles located in inner rows */
  update_position_inner_rows(eta);

  /* unpack data to row = 1*/
  MPI_Wait(&reqs[0], &stats[0]);
  int size_from_pre;
  MPI_Get_count(&stats[0], MPI_DOUBLE, &size_from_pre);
  unpack(1, buff_from_pre, size_from_pre);

  /* unpack data to row = nrows-2 */
  MPI_Wait(&reqs[2], &stats[2]);
  int size_from_next;
  MPI_Get_count(&stats[2], MPI_DOUBLE, &size_from_next);
  unpack(nrows - 2, buff_from_next, size_from_next);

  /* wait until finishing sending */
  MPI_Wait(&reqs[1], &stats[1]);
  MPI_Wait(&reqs[3], &stats[3]);

  /* free memory */
  delete[] buff_from_next;
  delete[] buff_from_pre;
  delete[] buff_to_next;
  delete[] buff_to_pre;
}

void SubDomain::one_step_MPI(double eta, int t) {
  update_velocity_MPI();
  update_position_MPI(eta);
  //remove_ghost_particle(0);
  //remove_ghost_particle(nrows - 1);
  if (t % 100 == 0) {
  int n2 = count_valid_particle();
  cout << "rank" << myrank << "\tN = " << n2 << "\tt = " << t << endl;
  }
}

