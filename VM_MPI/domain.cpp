#include "domain.h"
#include <iostream>
#include <cmath>
using namespace std;

/****************************************************************************/
/*                   Basic class for subdomain                              */
/****************************************************************************/
BasicDomain::BasicDomain(double eta, double eps, double rho0,
                         double Lx0, double Ly0, unsigned long long seed):
                         Lx(Lx0), Ly(Ly0) {
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  pre_rank = myrank == 0 ? tot_rank - 1 : myrank - 1;
  next_rank = myrank == tot_rank - 1 ? 0 : myrank + 1;

  double dLy = Ly / tot_rank;
  yl = myrank * dLy - 1;
  yh = (myrank + 1) * dLy + 1;
  ncols = int(Lx);
  nrows = int(dLy) + 2;
  cell = NULL;

  particle = NULL;
  MAX_PAR = end_pos = MAX_BUF_SIZE = nPar_per_row = nPar_per_task = 0;

  myran = new Ran(seed + myrank);
  if (myrank == 0) {
    char fname[100];
    snprintf(fname, 100, "p_%g_%g_%g_%g_%g_%llu_n%d.dat",
             eta, eps, rho0, Lx, Ly, seed, tot_rank);
    fout_phi.open(fname);
  }

  cout << "subdomain " << myrank << "\tpre rank " << pre_rank;
  cout << "\tnext rank " << next_rank << endl;
  cout << "\tncols = " << ncols << "\tnrows = " << nrows << endl;
}

void BasicDomain::create_particle_random(int nPar, double multiple) {
  MAX_PAR = int(nPar * multiple);
  particle = new Node[MAX_PAR];
  end_pos = nPar;
  Node::ini_random(particle, nPar, myran, Lx, yl, yh);
  create_cell_list();
  nPar_per_row = nPar / (nrows - 2);
  nPar_per_task = nPar;
  MAX_BUF_SIZE = nPar_per_row * 10 * 4;
}

void BasicDomain::create_from_snap(const string & filename) {
  Node::ini_from_snap(&particle, end_pos, MAX_PAR, filename,
                      Lx, Ly, tot_rank, myrank);
  create_cell_list();
  nPar_per_row = end_pos / (nrows - 2);
  nPar_per_task = end_pos;
  MAX_BUF_SIZE = nPar_per_row * 10 * 4;
}

void BasicDomain::output(int t) {
  if (t % 100 == 0) {
    double sv[2];
    int my_num;
    Node::sum_v(particle, end_pos, my_num, sv[0], sv[1]);
    double tot_v[2];
    int *num = new int[tot_rank];
    MPI_Gather(&my_num, 1, MPI_INT, num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Reduce(sv, tot_v, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      int sum_num = 0;
      for (int i = 0; i < tot_rank; i++) sum_num += num[i];
      double phi = sqrt(tot_v[0] * tot_v[0] + tot_v[1] * tot_v[1]) / sum_num;
      double theta = atan2(tot_v[1], tot_v[0]);
      fout_phi << t << "\t" << phi << "\t" << theta << "\t" << sum_num;
      for (int i = 0; i < tot_rank; i++) {
        fout_phi << "\t" << num[i] - sum_num / tot_rank;
      }
      fout_phi << endl;
    }
    delete[] num;
  }
}

void BasicDomain::create_cell_list() {
  for (int i = 0; i < end_pos; i++) {
    int idx = particle[i].cell_idx(yl, ncols);
    cell[idx].push_front(&particle[i]);
  }
}

void BasicDomain::update_position_edge_row(int row, double eta) {
  if (row != 1 && row != nrows - 2) {
    cout << "Error, row = " << row << " should be equal to "
      << 1 << " or " << nrows - 2 << endl;
    exit(1);
  }
  std::vector<Node *> ghost;
  ghost.reserve(MAX_BUF_SIZE / 8);
  for (int col = 0; col < ncols; col++) {
    int i = col + row * ncols;
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        double noise = eta * 2 * PI * (myran->doub() - 0.5);
        curNode->update_coor(noise, Lx, yl);
        if (curNode->y < yl + 1 || curNode->y >= yh - 1) {
          curNode->is_ghost = true;
          ghost.push_back(curNode);
        } else {
          curNode->is_ghost = false;
        }
        curNode->new_arrival = false;
        curNode = curNode->next;
      } while (curNode);
    }
  }
  for (int i = 0, size = ghost.size(); i < size; i++) {
    int c_idx = ghost[i]->cell_idx(yl, ncols);
    cell[c_idx].push_front(ghost[i]);
  }
}

void BasicDomain::update_position_inner_rows(double eta) {
  for (int row = 1; row < nrows - 1; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      cell[col + j].head = NULL;
      cell[col + j].size = 0;
    }
  }
  for (int i = 0; i < end_pos; i++) {
    if (!particle[i].is_empty) {
      if (!particle[i].is_moved) {
        double noise = eta * 2.0 * PI * (myran->doub() - 0.5);
        particle[i].update_coor(noise, Lx, yl);
        particle[i].is_ghost = false;
      }
      if (!particle[i].is_ghost) {
        int idx = particle[i].cell_idx(yl, ncols);
        cell[idx].push_front(&particle[i]);
      }
      particle[i].is_moved = false;   /* set all occupied particles unmoved */
    }
  }
}

void BasicDomain::pack(int row, double * buf, int &buf_size) {
  double dy = 0;
  if (myrank == 0 && row < nrows / 2) {
    dy = Ly;
  } else if (myrank == tot_rank - 1 && row > nrows / 2) {
    dy = -Ly;
  }
  int MAX_SIZE = buf_size;
  int pos = 0;
  int j = row * ncols;
  for (int col = 0; col < ncols; col++) {
    int i = col + j;
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        if (!curNode->new_arrival) {
          buf[pos++] = curNode->x;
          buf[pos++] = curNode->y + dy;
          buf[pos++] = curNode->vx;
          buf[pos++] = curNode->vy;
        } else {
          curNode->new_arrival = false;
        }
        curNode = curNode->next;
        if (pos >= MAX_SIZE) {
          cout << "Error, pos = " << pos << " is larger than MAX_BUFF_SIZE = "
            << MAX_SIZE << endl;
          exit(1);
        }
      } while (curNode);
    }
  }
  buf_size = pos;
}

void BasicDomain::unpack(int row, const double *buf, int buf_size) {
  bool ghost;
  if (row == 0 || row == nrows - 1) {
    ghost = true;
  } else if (row > 0 && row < nrows - 1) {
    ghost = false;
  } else {
    cout << "Failed to unpack, row = " << row << endl;
    exit(1);
  }
  int nPar = buf_size / 4;
  for (int i = 0, j = row * ncols; i < nPar; i++) {
    double x = buf[i * 4];
    double y = buf[i * 4 + 1];
    double vx = buf[i * 4 + 2];
    double vy = buf[i * 4 + 3];
    int cell_idx = int(x) + j;
    if (int(y - yl) != row) {
      cout << "Error, y = " << y << "\tint(y - yl) = "
        << int(y - yl) << "\trow = " << row << endl;
      exit(1);
    }
    Node *p = NULL;
    if (!empty_pos.empty()) {
      p = empty_pos.top();
      empty_pos.pop();
    } else {
      p = particle + end_pos;
      end_pos++;
    }
    p->x = x;
    p->y = y;
    p->vx = p->vx0 = vx;
    p->vy = p->vy0 = vy;
    p->is_moved = false;
    p->is_empty = false;
    p->is_ghost = ghost;
    p->new_arrival = true;
    cell[cell_idx].push_front(p);
  }
}

void BasicDomain::update_position(double eta) {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];

  for (int i = 0; i < 4; i++) {
    buf[i] = new double[MAX_BUF_SIZE];
    buf_size[i] = MAX_BUF_SIZE;
  }

  /* pre_rank -> myrank */
  recv(buf[0], buf_size[0], pre_rank, 12, &req[0]);
  /* myrank -> next_rank */
  update_position_edge_row(nrows - 2, eta);
  send(nrows - 1, buf[1], buf_size[1], next_rank, 12, &req[1]);
  /* myrank <- next_rank */
  recv(buf[2], buf_size[2], next_rank, 21, &req[2]);
  /* pre_rank <- myrank */
  update_position_edge_row(1, eta);
  send(0, buf[3], buf_size[3], pre_rank, 21, &req[3]);

  /* update position of inner rows */
  update_position_inner_rows(eta);

  /* pre_rank -> myrank */
  accept(1, buf[0], &buf_size[0], &req[0], &stat[0]);

  /* myrank <- next_rank */
  accept(nrows - 2, buf[2], &buf_size[2], &req[2], &stat[2]);

  MPI_Wait(&req[1], &stat[1]);
  MPI_Wait(&req[3], &stat[3]);
  for (int i = 0; i < 4; i++) {
    delete[] buf[i];
  }
}

/****************************************************************************/
/*                   Static domain decomposition                            */
/****************************************************************************/
StaticDomain::StaticDomain(double eta, double eps, double rho0,
                           double Lx0, double Ly0, unsigned long long seed):
                           BasicDomain(eta, eps, rho0, Lx0, Ly0, seed) {
  cell = new Cell[ncols * nrows];
  Cell::find_all_neighbor(cell, ncols, nrows);
}

StaticDomain::~StaticDomain() {
  delete[] cell;
  delete[] particle;
  fout_phi.close();
}

void StaticDomain::update_velocity() {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];
  for (int i = 0; i < 4; i++) {
    buf[i] = new double[MAX_BUF_SIZE];
    buf_size[i] = MAX_BUF_SIZE;
  }

  /* pre_rank -> myrank */
  MPI_Irecv(buf[0], buf_size[0], MPI_DOUBLE,
            pre_rank, 13, MPI_COMM_WORLD, &req[0]);
  /* myrank -> next_rank */
  send(nrows - 2, buf[1], buf_size[1], next_rank, 13, &req[1]);
  /* myrank <- next_rank */
  MPI_Irecv(buf[2], buf_size[2], MPI_DOUBLE,
            next_rank, 31, MPI_COMM_WORLD, &req[2]);
  /* pre_rank <- myrank */
  send(1, buf[3], buf_size[3], pre_rank, 31, &req[3]);

  for (int row = 1; row < nrows - 2; row++) {
    Cell::update_velocity_inner_row(cell + row * ncols, ncols, Lx);
  }

  /* pre_rank -> myrank */
  accept(0, buf[0], &buf_size[0], &req[0], &stat[0]);
  Cell::update_velocity_by_row(0, cell, ncols, Lx);
  Cell::clear_row(0, ncols, cell, empty_pos);
  /* myrank <- next_rank */
  accept(nrows - 1, buf[2], &buf_size[2], &req[2], &stat[2]);
  Cell::update_velocity_by_row(nrows - 2, cell, ncols, Lx);
  Cell::clear_row(nrows - 1, ncols, cell, empty_pos);

  MPI_Wait(&req[1], &stat[1]);
  MPI_Wait(&req[3], &stat[3]);

  for (int i = 0; i < 4; i++) {
    delete[] buf[i];
  }
}

void StaticDomain::one_step(double eta, int t) {
  update_velocity();
  update_position(eta);
  output(t);
}

/****************************************************************************/
/*                   Dynamic domain decomposition                           */
/****************************************************************************/

DynamicDomain::DynamicDomain(double eta, double eps, double rho0,
                             double Lx0, double Ly0, unsigned long long seed):
                             BasicDomain(eta, eps, rho0, Lx0, Ly0, seed) {
  CELL_BUF_SIZE = ncols * (nrows - 2) * 3;
  cell_buf = new Cell[CELL_BUF_SIZE];
  int first_row = nrows - 2;
  cell = cell_buf + first_row * ncols;
  Cell::find_all_neighbor(cell, ncols, nrows);
  row_offset[0] = row_offset[1] = 0;
}

DynamicDomain::~DynamicDomain() {
  delete[] cell_buf;
  delete[] particle;
  fout_phi.close();
}

void DynamicDomain::rearrange_domain(int t) {
  int my_nPar[3];
  Cell::get_nPar(cell, ncols, nrows, my_nPar);
  int *recvbuf = new int[tot_rank * 3];
  int *sendbuf = new int[tot_rank * 2];

  MPI_Gather(&my_nPar, 3, MPI_INT, recvbuf, 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    for (int i = 0; i < tot_rank * 2; i++)
      sendbuf[i] = 0;
    int threshold = 2 * nPar_per_row;
    for (int i = 0; i < tot_rank - 1; i++) {
      int j = 3 * i + 1;
      int k = 2 * i + 1;
      int dn = recvbuf[j] - nPar_per_task;
      if (dn > threshold) {
        if (recvbuf[j] > recvbuf[j + 3] || dn > 2 * threshold) {
          recvbuf[j] -= recvbuf[j + 1];
          recvbuf[j + 3] += recvbuf[j + 1];
          sendbuf[k] = -1;
          sendbuf[k + 1] = 1;
        }
      } else if (dn < -threshold) {
        if (recvbuf[j] < recvbuf[j + 3] || dn < -2 * threshold) {
          recvbuf[j] += recvbuf[j + 1];
          recvbuf[j + 3] -= recvbuf[j + 1];
          sendbuf[k] = 1;
          sendbuf[k + 1] = -1;
        }
      }
    }
  }
  MPI_Scatter(sendbuf, 2, MPI_INT, row_offset, 2, MPI_INT, 0, MPI_COMM_WORLD);
  delete[] recvbuf;
  delete[] sendbuf;
}

void DynamicDomain::update_velocity() {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];


  for (int i = 0; i < 4; i++) {
    buf[i] = new double[MAX_BUF_SIZE];
    buf_size[i] = MAX_BUF_SIZE;
  }
  /* pack and transfer the data */
  {
    int i = 0;
    /* pre_rank -> myrank */
    for (int j = 0; j <= row_offset[0]; j++) {
      recv(buf[i], buf_size[i], pre_rank, 12 + j, &req[i]);
      i++;
    }
    /* myrank -> next_rank */
    for (int j = 0; j >= row_offset[1]; j--) {
      send(nrows - 2 + j, buf[i], buf_size[i], next_rank, 12 - j, &req[i]);
      i++;
    }
    /* myrank <- next_rank */
    for (int j = 0; j <= row_offset[1]; j++) {
      recv(buf[i], buf_size[i], next_rank, 21 + j * 10, &req[i]);
      i++;
    }
    /* pre_rank <- myrank */
    for (int j = 0; j >= row_offset[0]; j--) {
      send(1 - j, buf[i], buf_size[i], pre_rank, 21 - j * 10, &req[i]);
      i++;
    }
  }
  /* update velocity for inner rows */
  for (int row = 1; row < nrows - 2; row++)
    Cell::update_velocity_inner_row(cell + row * ncols, ncols, Lx);

  /* update the size of the cell list */
  Cell::resize(&cell, ncols, nrows, row_offset, empty_pos);
  yl -= row_offset[0];
  yh += row_offset[1];

  /* unpack the received data and update velocity for remaining rows */
  {
    int i = 0;
    for (int j = 0; j <= row_offset[0]; j++) {
      int row = row_offset[0] == 0 ? 0 : 1 - j;
      accept(row, buf[i], &buf_size[i], &req[i], &stat[i]);
      Cell::update_velocity_by_row(row, cell, ncols, Lx);
      i++;
    }
    Cell::clear_row(0, ncols, cell, empty_pos);

    for (int j = 0; j >= row_offset[1]; j--) i++;

    for (int j = 0; j <= row_offset[1]; j++) {
      int row = row_offset[1] == 0 ? nrows - 1 : nrows - 2 + j;
      accept(row, buf[i], &buf_size[i], &req[i], &stat[i]);
      Cell::update_velocity_by_row(row - 1, cell, ncols, Lx);
      i++;
    }
    Cell::clear_row(nrows - 1, ncols, cell, empty_pos);
  }

  MPI_Waitall(4, req, stat);
  for (int i = 0; i < 4; i++) {
    delete[] buf[i];
    buf[i] = NULL;
  }
}

void DynamicDomain::one_step(double eta, int t) {
  rearrange_domain(t);
  update_velocity();
  update_position(eta);
  output(t);
}
