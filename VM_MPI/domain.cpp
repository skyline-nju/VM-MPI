#include "domain.h"
#include <iostream>
#include <cmath>
using namespace std;

/****************************************************************************/
/*                        Common functions                                  */
/****************************************************************************/

void remove_ghost_particle(int row, int ncols, Cell *cell, 
                           stack<Node *> &empty_pos) {
  for (int col = 0, j = row * ncols; col < ncols; col++) {
    int idx = col + j;
    if (cell[idx].head) {
      Node *curNode = cell[idx].head;
      do {
        curNode->is_empty = true;
        empty_pos.push(curNode);
        curNode = curNode->next;
      } while (curNode);
      cell[idx].head = NULL;
      cell[idx].size = 0;
    }
  }
}

void sum_velocity(const Node *par, int end_pos,
                  int &npar, double &svx, double &svy) {
  svx = 0;
  svy = 0;
  npar = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!par[i].is_empty && !par[i].is_ghost) {
      svx += par[i].vx;
      svy += par[i].vy;
      npar++;
    }
  }
}

int get_count_par_num(const Cell *cell, int ncols, int nrows) {
  int count = 0;
  for (int row = 1; row < nrows - 1; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      count += cell[col + j].size;
    }
  }
  return count;
}

void get_count_par_num(const Cell *cell, int ncols, int nrows, int *res) {
  res[0] = res[1] = res[2] = 0;
  for (int col = 0; col < ncols; col++) {
    res[0] += cell[col + ncols].size;
  }
  for (int row = 2; row < nrows - 2; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      res[1] += cell[col + j].size;
    }
  }
  for (int col = 0, j = (nrows - 2) * ncols; col < ncols; col++) {
    res[2] += cell[col + j].size;
  }
  res[1] += (res[0] + res[2]);
}

void create_cell_list(Cell *cell, int ncols, double yl,
                      Node *par, int nPar) {
  for (int i = 0; i < nPar; i++) {
    int idx = par[i].cell_idx(yl, ncols);
    cell[idx].push_front(&par[i]);
  }
}

void update_position_inner_rows(double eta, Ran *myran,
                                Node *par, int end_pos,
                                Cell *cell, int ncols, int nrows, 
                                double Lx, double yl) {
  for (int row = 1; row < nrows - 1; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      cell[col + j].head = NULL;
      cell[col + j].size = 0;
    }
  }
  for (int i = 0; i < end_pos; i++) {
    if (!par[i].is_empty) {
      if (!par[i].is_moved) {
        double noise = eta * 2.0 * PI * (myran->doub() - 0.5);
        par[i].update_coor(noise, Lx, yl);
      }
      par[i].is_moved = false;
      if (!par[i].is_ghost) {
        int idx = par[i].cell_idx(yl, ncols);
        cell[idx].push_front(&par[i]);
      }
    }
  }
}

void update_position_edge_row(int row, double eta, Ran *myran,
                              Cell *cell,int ncols, int nrows,
                              double Lx, double yl, double yh,
                              int MAX_BUFF_SIZE) {
  if (row != 1 && row != nrows - 2) {
    cout << "Error, row = " << row << " should be equal to "
         << 1 << " or " << nrows - 2 << endl;
    exit(1);
  }
  std::vector<Node *> ghost;
  ghost.reserve(MAX_BUFF_SIZE / 8);
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
        }
        curNode = curNode->next;
      } while (curNode);
    }
  }
  for (int i = 0, size = ghost.size(); i < size; i++) {
    int c_idx = ghost[i]->cell_idx(yl, ncols);
    cell[c_idx].push_front(ghost[i]);
  }
}

void pack(int row, double *buff, int &buff_size,
          Cell *cell, int nrows, int ncols, double Ly,
          int myrank, int tot_rank, int MAX_BUFF_SIZE) {
  double dy = 0;
  if (myrank == 0 && row < nrows / 2) {
    dy = Ly;
  } else if (myrank == tot_rank - 1 && row > nrows / 2) {
    dy = -Ly;
  }
  int pos = 0;
  int j = row * ncols;
  for (int col = 0; col < ncols; col++) {
    int i = col + j;
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        if (!curNode->new_arrival) {
          buff[pos++] = curNode->x;
          buff[pos++] = curNode->y + dy;
          buff[pos++] = curNode->vx;
          buff[pos++] = curNode->vy;
        } else {
          curNode->new_arrival = false;
        }
        curNode = curNode->next;
        if (pos >= MAX_BUFF_SIZE) {
          cout << "Error, pos = " << pos << " is larger than MAX_BUFF_SIZE = "
               << MAX_BUFF_SIZE << endl;
          exit(1);
        }
      } while (curNode);
    }
  }
  buff_size = pos;
}

void unpack(int row, const double *buff, int buff_size,
            Cell *cell, int nrows, int ncols, double yl,
            Node *par, int &end_pos, stack<Node *> &empty_pos) {
  bool ghost;
  if (row == 0 || row == nrows - 1) {
    ghost = true;
  } else if (row > 0 && row < nrows - 1) {
    ghost = false;
  } else {
    cout << "Failed to unpack, row = " << row << endl;
    exit(1);
  }
  int nPar = buff_size / 4;
  for (int i = 0, j = row * ncols; i < nPar; i++) {
    double x = buff[i * 4];
    double y = buff[i * 4 + 1];
    double vx = buff[i * 4 + 2];
    double vy = buff[i * 4 + 3];
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
      p = par + end_pos;
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

/****************************************************************************/
/*                   Static domain decomposition                            */
/****************************************************************************/

StaticDomain::StaticDomain(double Lx_domain, double Ly_domain, int ntask,
                           int rank, unsigned long long seed, double eta,
                           double eps, double rho0){
  Lx = Lx_domain;
  Ly = Ly_domain;
  double dLy = Ly / ntask;
  yl = rank * dLy - 1;
  yh = (rank + 1) * dLy + 1;
  ncols = int(Lx);
  nrows = int(dLy) + 2;
  int tot_cells = ncols * nrows;
  cell = new Cell[tot_cells];
  Cell::find_all_neighbor(cell, ncols, nrows);
  myran = new Ran(seed + rank);
  tot_rank = ntask;
  myrank = rank;
  pre_rank = rank == 0 ? ntask - 1 : rank - 1;
  next_rank = rank == ntask - 1 ? 0 : rank + 1;
  cout << "subdomain " << rank << "\tpre rank " << pre_rank;
  cout << "\tnext rank " << next_rank << endl;
  cout << "\tncols = " << ncols << "\tnrows = " << nrows << endl;

  if (myrank == 0) {
    char fname[100];
    snprintf(fname, 100, "p_%g_%g_%g_%g_%g_%llu_n%d.dat",
             eta, eps, rho0, Lx, Ly, seed, tot_rank);
    fout_phi.open(fname);
  }
}

StaticDomain::~StaticDomain() {
  delete[] cell;
  delete[] particle;
  fout_phi.close();
}

void StaticDomain::create_particle_random(int nPar) {
  MAX_PAR = nPar * 2;
  particle = new Node[MAX_PAR];
  Node::ini_random(particle, nPar, myran, Lx, yl, yh);
  create_cell_list(cell, ncols, yl, particle, end_pos);
  MAX_BUFF_SIZE = nPar / (nrows - 2) * 10 * 4;
}

void StaticDomain::create_from_snap(const string filename) {
  Node::ini_from_snap(&particle, end_pos, MAX_PAR, filename,
                      Lx, Ly, tot_rank, myrank);
  create_cell_list(cell, ncols, yl, particle, end_pos);
  cout << "rank = " << myrank << "\tN = " << end_pos << endl;
  MAX_BUFF_SIZE = end_pos / (nrows - 2) * 10 * 4;
}

void StaticDomain::comm_start(int source_row, int &dest_row,
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
    cout << "Wrong row when starting communication\n";
    exit(1);
  }
  double *send_buff = new double[MAX_BUFF_SIZE];
  double *recv_buff = new double[MAX_BUFF_SIZE];
  int send_buff_size;
  pack(source_row, send_buff, send_buff_size, cell, nrows, ncols, Ly,
       myrank, tot_rank, MAX_BUFF_SIZE);
  MPI_Irecv(
    recv_buff, MAX_BUFF_SIZE, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, rreq);
  MPI_Isend(
    send_buff, send_buff_size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, sreq);
  *sbuff = send_buff;
  *rbuff = recv_buff;
}

void StaticDomain::comm_end(int dest_row, double *sbuff, double *rbuff,
                            MPI_Request *sreq, MPI_Request *rreq) {
  MPI_Status stat;
  MPI_Wait(rreq, &stat);
  int recv_buff_size;
  MPI_Get_count(&stat, MPI_DOUBLE, &recv_buff_size);
  unpack(dest_row, rbuff, recv_buff_size, cell, nrows, ncols, yl,
         particle, end_pos, empty_pos);
  delete[] rbuff;
  MPI_Wait(sreq, &stat);
  delete[] sbuff;
}

void StaticDomain::update_velocity_MPI() {
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

  for (int row = 1; row < nrows - 2; row++) {
    Cell::update_velocity_inner_row(cell + row * ncols, ncols, Lx);
  }

  for (int i = 0; i < sreq.size(); i++) {
    comm_end(dest_row[i], sbuff[i], rbuff[i], &sreq[i], &rreq[i]);
    int row = dest_row[i] < nrows / 2 ? dest_row[i] : dest_row[i] - 1;
    if (row == 0) {
      Cell::update_velocity_bottom_row(cell, ncols, Lx);
    }
    else{
      Cell::update_velocity_inner_row(cell + ncols * row, ncols, Lx);
    }
    remove_ghost_particle(dest_row[i], ncols, cell, empty_pos);
  }
}

void StaticDomain::update_position_MPI(double eta) {
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
    update_position_edge_row(row, eta, myran, cell, ncols, nrows,
                             Lx, yl, yh, MAX_BUFF_SIZE);
    comm_start(src_row[i], dest_row[i], &sbuff[i], &rbuff[i],
               &sreq[i], &rreq[i]);
  }
  update_position_inner_rows(eta, myran, particle, end_pos,
                             cell, ncols, nrows, Lx, yl);
  for (int i = 0; i < sreq.size(); i++) {
    comm_end(dest_row[i], sbuff[i], rbuff[i], &sreq[i], &rreq[i]);
  }
}

void StaticDomain::one_step_MPI(double eta, int t) {
  update_velocity_MPI();
  update_position_MPI(eta);
  output(t);
}

void StaticDomain::output(int t) {
  if (t % 100 == 0) {
    double sv[2];
    int sub_N;
    sum_velocity(particle, end_pos, sub_N, sv[0], sv[1]);
    double tot_v[2];
    int *num = new int[tot_rank];
    MPI_Gather(&sub_N, 1, MPI_INT, num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Reduce(sv, tot_v, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      int totN = 0;
      double sum_N2 = 0;
      for (int i = 0; i < tot_rank; i++) {
        totN += num[i];
        sum_N2 += num[i] * num[i];
      }
      double mean_N = totN / tot_rank;
      double phi = sqrt(tot_v[0] * tot_v[0] + tot_v[1] * tot_v[1]) / totN;
      double theta = atan2(tot_v[1], tot_v[0]);
      fout_phi << t << "\t" << phi << "\t" << theta << "\t" << totN;
      for (int i = 0; i < tot_rank; i++) {
        fout_phi << "\t" << num[i] - totN / tot_rank;
      }
      fout_phi << endl;
    }
    delete[] num;
  }
}

/****************************************************************************/
/*                   Dynamic domain decomposition                           */
/****************************************************************************/

DynamicDomain::DynamicDomain(double Lx_domain, double Ly_domain, int ntask,
                             int rank, unsigned long long seed, double eta,
                             double eps, double rho0) {
  Lx = Lx_domain;
  Ly = Ly_domain;
  double dLy = Ly / ntask;
  yl = rank * dLy - 1;
  yh = (rank + 1) * dLy + 1;
  ncols = int(Lx);
  nrows = int(dLy) + 2;
  int first_row = nrows - 2;
  MAX_CELL = ncols * (nrows -2) * 3;
  cell_buff = new Cell[MAX_CELL];
  cell = cell_buff + first_row * ncols;
  Cell::find_all_neighbor(cell, ncols, nrows);
  row_offset[0] = row_offset[1] = 0;
  particle = NULL;
  myran = new Ran(seed + rank);
  tot_rank = ntask;
  myrank = rank;
  pre_rank = rank == 0 ? ntask - 1 : rank - 1;
  next_rank = rank == ntask - 1 ? 0 : rank + 1;
  cout << "subdomain " << rank << "\tpre rank " << pre_rank;
  cout << "\tnext rank " << next_rank << endl;
  cout << "\tncols = " << ncols << "\tnrows = " << nrows << endl;

  if (myrank == 0) {
    char fname[100];
    snprintf(fname, 100, "p_%g_%g_%g_%g_%g_%llu_n%d.dat",
      eta, eps, rho0, Lx, Ly, seed, tot_rank);
    fout_phi.open(fname);
  }
}

DynamicDomain::~DynamicDomain() {
  delete[] cell_buff;
  delete[] particle;
  fout_phi.close();
}

void DynamicDomain::create_particle_random(int nPar) {
  MAX_PAR_BUFF = nPar * 2;
  particle = new Node[MAX_PAR_BUFF];
  end_pos = nPar;
  Node::ini_random(particle, nPar, myran, Lx, yl, yh);
  create_cell_list(cell, ncols, yl, particle, end_pos);
  nPar_per_row = nPar / (nrows - 2);
  nPar_per_task = nPar;
  MAX_BUFF_SIZE = nPar_per_row * 10 * 4;
}

void DynamicDomain::create_from_snap(const std::string filename) {
  Node::ini_from_snap(&particle, end_pos, MAX_PAR_BUFF, filename,
                      Lx, Ly, tot_rank, myrank);
  create_cell_list(cell, ncols, yl, particle, end_pos);
  nPar_per_row = end_pos / (nrows - 2);
  nPar_per_task = end_pos;
  MAX_BUFF_SIZE = nPar_per_row * 10 * 4;
}

void DynamicDomain::rearrange() {
  int my_nPar[3];
  get_count_par_num(cell, ncols, nrows, my_nPar);
  int *recvbuf = new int[tot_rank];
  int *sendbuf = new int[tot_rank * 2];

  MPI_Gather(&my_nPar, 3, MPI_INT, recvbuf, 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    sendbuf[0] = 0;
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
        } else {
          sendbuf[k] = 0;
          sendbuf[k + 1] = 0;
        }
      } else if (dn < -threshold) {
        if (recvbuf[j] < recvbuf[j + 3] || dn < -2 * threshold) {
          recvbuf[j] += recvbuf[j + 1];
          recvbuf[j + 3] -= recvbuf[j + 1];
          sendbuf[k] = 1;
          sendbuf[k + 1] = -1;
        } else {
          sendbuf[k] = 0;
          sendbuf[k + 1] = 0;
        }
      }
    }
    sendbuf[2 * tot_rank - 1] = 0;
  }
  MPI_Scatter(sendbuf, 2, MPI_INT, row_offset, 2, MPI_INT, 0, MPI_COMM_WORLD);
  delete[] recvbuf;
  delete[] sendbuf;
}

void DynamicDomain::update_velocity() {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buff[4];
  int buff_size[4];
  int src_row[4];
  int dest_row[4];

  for (int i = 0; i < 4; i++) {
    buff[i] = new double[MAX_BUFF_SIZE];
    buff_size[i] = MAX_BUFF_SIZE;
  }
  
  /* pack and transfer the data */
  {
    int i = 0;
    /* pre_rank -> myrank -> next_rank */
    if (row_offset[0] == 1) {
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                pre_rank, 12, MPI_COMM_WORLD, &req[i]);
      i++;
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                pre_rank, 13, MPI_COMM_WORLD, &req[i]);
      i++;
    } else if (row_offset[0] == 0) {
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                pre_rank, 12, MPI_COMM_WORLD, &req[i]);
      i++;
    }
    if (row_offset[1] == -1) {
      pack(nrows - 2, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                next_rank, 12, MPI_COMM_WORLD, &req[i]);
      i++;
      pack(nrows - 3, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                next_rank, 13, MPI_COMM_WORLD, &req[i]);
      i++;
    } else if (row_offset[1] == 0) {
      pack(nrows - 2, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                next_rank, 12, MPI_COMM_WORLD, &req[i]);
      i++;
    }
    /* pre_rank <- myrank <- next_rank */
    if (row_offset[1] == 1) {
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                next_rank, 21, MPI_COMM_WORLD, &req[i]);
      i++;
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                next_rank, 31, MPI_COMM_WORLD, &req[i]);
      i++;
    } else if (row_offset[1] == 0) {
      MPI_Irecv(buff[i], MAX_BUFF_SIZE, MPI_DOUBLE,
                next_rank, 21, MPI_COMM_WORLD, &req[i]);
      i++;
    }
    if (row_offset[0] == -1) {
      pack(1, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                pre_rank, 21, MPI_COMM_WORLD, &req[i]);
      i++;
      pack(2, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                pre_rank, 31, MPI_COMM_WORLD, &req[i]);
      i++;
    } else if (row_offset[1] == 0) {
      pack(1, buff[i], buff_size[i],
           cell, nrows, ncols, Ly, myrank, tot_rank, MAX_BUFF_SIZE);
      MPI_Isend(buff[i], buff_size[i], MPI_DOUBLE,
                pre_rank, 21, MPI_COMM_WORLD, &req[i]);
      i++;
    }
  }
  
  /* update velocity for inner rows */
  {
    int beg_row = 1;
    int end_row = nrows - 2;
    for (int row = beg_row; row < end_row; row++)
      Cell::update_velocity_inner_row(cell + row * ncols, ncols, Lx);
  }

  /* unpack the data and update velocity for remaining rows */
  {
    /* update the cell list */
    Cell::update_neighbor(&cell, ncols, nrows, row_offset);
    yl -= row_offset[0];
    yh += row_offset[1];

  }

  for (int i = 0; i < 4; i++) {
    delete[] buff[i];
  }
}

void DynamicDomain::update_position(double eta) {
}
