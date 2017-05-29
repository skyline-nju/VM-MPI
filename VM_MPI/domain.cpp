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
      cell[col + j].head = NULL;
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
      cell[idx].head = NULL;
      cell[idx].size = 0;
    }
  }
}


int SubDomain::count_valid_particle() const {
  int count = 0;
  for (int i = 0, size = particle.size(); i < size; i++) {
    if (!particle[i].is_empty && !particle[i].is_ghost)
      count += 1;
  }
  return count;
}

void SubDomain::sum_velocity(double &svx, double &svy, int &npar) const {
  svx = 0;
  svy = 0;
  npar = 0;
  for (int i = 0, size = particle.size(); i < size; i++) {
    if (!particle[i].is_empty && !particle[i].is_ghost) {
      svx += particle[i].vx;
      svy += particle[i].vy;
      npar++;
    }
  }
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
  for (int i = 0; i < nPar; i++) {
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

void SubDomain::comm_start(int source_row, int &dest_row,
                           double **sbuff, double **rbuff,
                           MPI_Request *sreq, MPI_Request *rreq) {
  int source, dest, tag = source_row;
  if (source_row == nrows - 2 || source_row == nrows - 1) {
    source = pre_rank;
    dest = next_rank;
    dest_row = source_row - nrows + 2;
  } else if (source_row == 0 || source_row == 1) {
    source = next_rank;
    dest = pre_rank;
    dest_row = source_row + nrows - 2;
  } else {
    cout << "Wrong row when start communication\n";
    exit(1);
  }
  double *send_buff = new double[MAX_BUFF_SIZE];
  double *recv_buff = new double[MAX_BUFF_SIZE];
  int send_buff_size;
  pack(source_row, send_buff, send_buff_size);
  MPI_Irecv(
    recv_buff, MAX_BUFF_SIZE, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, rreq);
  MPI_Isend(
    send_buff, send_buff_size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, sreq);
  *sbuff = send_buff;
  *rbuff = recv_buff;
}

void SubDomain::comm_end(int dest_row, double *sbuff, double *rbuff,
                         MPI_Request *sreq, MPI_Request *rreq) {
  MPI_Status stat;
  MPI_Wait(rreq, &stat);
  int recv_buff_size;
  MPI_Get_count(&stat, MPI_DOUBLE, &recv_buff_size);
  unpack(dest_row, rbuff, recv_buff_size);
  delete[] rbuff;
  MPI_Wait(sreq, &stat);
  delete[] sbuff;
}

void SubDomain::update_velocity_MPI() {
  vector<MPI_Request> sreq(2);
  vector<MPI_Request> rreq(2);
  vector<double *> sbuff(2);
  vector<double *> rbuff(2);
  vector<int> dest_row(2);
  vector<int> source_row;
  source_row.push_back(nrows - 2);
  source_row.push_back(1);

  for (int i = 0; i < sreq.size(); i++) {
    comm_start(source_row[i], dest_row[i], &sbuff[i], &rbuff[i],
               &sreq[i], &rreq[i]);
  }

  update_velocity_inner_rows();

  for (int i = 0; i < sreq.size(); i++) {
    comm_end(dest_row[i], sbuff[i], rbuff[i], &sreq[i], &rreq[i]);
    int row = dest_row[i] < nrows / 2 ? dest_row[i] : dest_row[i] - 1;
    update_velocity_by_row(row);
    remove_ghost_particle(dest_row[i]);
  }
}

void SubDomain::update_position_MPI(double eta) {
  vector<MPI_Request> sreq(2);
  vector<MPI_Request> rreq(2);
  vector<double *> sbuff(2);
  vector<double *> rbuff(2);
  vector<int> dest_row(2);
  vector<int> src_row;
  src_row.push_back(nrows - 1);
  src_row.push_back(0);

  for (int i = 0; i < sreq.size(); i++) {
    int row = src_row[i] < nrows / 2 ? src_row[i] + 1 : src_row[i] - 1;
    update_position_edge_row(eta, row);
    comm_start(src_row[i], dest_row[i], &sbuff[i], &rbuff[i],
               &sreq[i], &rreq[i]);
  }

  update_position_inner_rows(eta);

  for (int i = 0; i < sreq.size(); i++) {
    comm_end(dest_row[i], sbuff[i], rbuff[i], &sreq[i], &rreq[i]);
  }
}

void SubDomain::one_step_MPI(double eta, int t, double &sum_phi, int &count) {
  update_velocity_MPI();
  update_position_MPI(eta);
  if (t % 100 == 0) {
    double sv[2];
    int N;
    sum_velocity(sv[0], sv[1], N);
    double tot_v[2];
    int *num = new int[tot_rank];
    MPI_Gather(&N, 1, MPI_INT, num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Reduce(sv, tot_v, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      int totN = 0;
      for (int i = 0; i < tot_rank; i++) {
        totN += num[i];
      }
      double phi = sqrt(tot_v[0] * tot_v[0] + tot_v[1] * tot_v[1]) / totN;
      if (t % 1000 == 0) {
        cout << "t = " << t << "\t";
        for (int i = 0; i < tot_rank; i++) {
          if (i == tot_rank - 1)
            cout << num[i] << " = " << totN;
          else
            cout << num[i] << " + ";
        }
        cout << "\tphi = " << phi;
        if (t >= 100000) {
          count++;
          sum_phi += phi;
          cout << "\t" << sum_phi / count << endl;
        } else {
          cout << endl;
        }
      }
    }
    delete[] num;
  }
}

