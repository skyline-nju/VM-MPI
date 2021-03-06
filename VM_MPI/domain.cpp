#include "domain.h"
#include <iostream>
#include <cmath>
using namespace std;

/****************************************************************************/
/*                   Basic class for subdomain                              */
/****************************************************************************/
BasicDomain::BasicDomain(const cmdline::parser &cmd) {
  eta = cmd.get<double>("eta");
  eps = cmd.get<double>("eps");
  rho0 = cmd.get<double>("rho0");
  Lx = cmd.get<double>("Lx");
  Ly = cmd.get<double>("Ly");
  tot_steps = cmd.get<int>("nstep");
  unsigned long long seed = cmd.get<unsigned long long>("seed");

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
  my_nPar = 0;
  MAX_PAR = end_pos = MAX_BUF_SIZE = nPar_per_row = nPar_per_task = 0;

  myran = new Ran(seed + myrank);
  disorder_free = eps > 0 ? false : true;
  cgs = NULL;
}

void BasicDomain::create_particle_random(const cmdline::parser &cmd,
                                         int tot_nPar, double magnification) {
  my_nPar = tot_nPar / tot_rank;
  MAX_PAR = int(my_nPar * magnification);
  particle = new Node[MAX_PAR];
  end_pos = my_nPar;
  Node::ini_random(particle, my_nPar, myran, Lx, yl, yh);
  create_cell_list();
  nPar_per_row = my_nPar / (nrows - 2);
  nPar_per_task = my_nPar;
  MAX_BUF_SIZE = nPar_per_row * 10 * 4;

  if (cmd.get<double>("cg_exp") > 0 || cmd.get<int>("cg_dt") > 0) {
    cout << "coarse grain" << endl;
    if (myrank == 0)
      mkdir("coarse");
    MPI_Barrier(MPI_COMM_WORLD);
    cgs = new CoarseGrainSnap(cmd, my_nPar * tot_rank);
  }
  if (myrank == 1) {
    fout << "create particles with random positions and velocities." << endl;
    if (cgs) {
      fout << "output coarse-grained snapshots with ";
      if (cmd.get<int>("cg_dt") > 0)
        fout << "dt = " << cmd.get<int>("cg_dt");
      if (cmd.get<double>("cg_exp") > 0)
        fout << ", exponent = " << cmd.get<double>("cg_exp");
      fout << endl;
    }
    fout << "-------- Run --------\n";
    fout << "time step\telapsed time\n";
  }
}

void BasicDomain::create_from_snap(const cmdline::parser & cmd,
                                   double magnification) {
  string fname = cmd.get<string>("fin");
  Node::ini_from_snap(&particle, end_pos, MAX_PAR, magnification, fname,
                      Lx, Ly, tot_rank, myrank);
  if (cmd.exist("tilt")) {
    double k = cmd.get<double>("tilt");
    for (int i = 0; i < end_pos; i++) {
      particle[i].x += k * particle[i].y;
      particle[i].x = fmod(particle[i].x, Lx);
    }
  }
  create_cell_list();
  nPar_per_row = end_pos / (nrows - 2);
  my_nPar = nPar_per_task = end_pos;
  MAX_BUF_SIZE = nPar_per_row * 10 * 4;

  if (cmd.get<double>("cg_exp") > 0 || cmd.get<int>("cg_dt") > 0) {
    cout << "coarse grain" << endl;
    if (myrank == 0)
      mkdir("coarse");
    MPI_Barrier(MPI_COMM_WORLD);
    cgs = new CoarseGrainSnap(cmd, my_nPar * tot_rank);
  }

  if (myrank == 1) {
    fout << "create particles by copying the snapshot: " << fname << endl;
    if (cgs) {
      fout << "output coarse-grained snapshots with ";
      if (cmd.get<int>("cg_dt") > 0)
        fout << "dt = " << cmd.get<int>("cg_dt");
      if (cmd.get<double>("cg_exp") > 0)
        fout << ", exponent = " << cmd.get<double>("cg_exp");
      fout << endl;
    }
    fout << "-------- Run --------\n";
    fout << "time step\telapsed time\n";
  }
}

void BasicDomain::output(int t) {
  if (cgs) {
    cgs->output(particle, end_pos, yl, nrows, t);
  }
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
      fout << t << "\t" << phi << "\t" << theta << "\t" << sum_num;
      for (int i = 0; i < tot_rank; i++) {
        fout << "\t" << num[i] - sum_num / tot_rank;
      }
      fout << endl;
    }
    delete[] num;
    if (t % 10000 == 0 && myrank == 1) {
      double dt = difftime(time(NULL), beg_time);
      int hour = int(dt / 3600);
      int min = int((dt - hour * 3600) / 60);
      double sec = dt - hour * 3600 - min * 60;
      fout << t << "\t" << hour << ":" << min << ":" << sec << endl;
    }
  }
}

void BasicDomain::ini_output(const string &tag, double eta, double eps,
                             double rho0 ,unsigned long long seed) {
  char fname[100];
  if (myrank == 0) {
    mkdir("phi");
    snprintf(fname, 100, "phi%sp_%g_%g_%g_%g_%g_%llu_n%d%s.dat",
             delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed, tot_rank, tag.c_str());
    fout.open(fname);
  } else if (myrank == 1) {
    mkdir("log");
    snprintf(fname, 100, "log%sl_%g_%g_%g_%g_%g_%llu_n%d%s.dat",
             delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed, tot_rank, tag.c_str());
    fout.open(fname);
    time(&beg_time);
    fout << "Started at " << asctime(localtime(&beg_time));
    fout << "-------- Parameters --------\n";
    fout << "eta: " << eta << endl;
    fout << "epsilon: " << eps << endl;
    fout << "density: " << rho0 << endl;
    fout << "Lx: " << Lx << endl;
    fout << "Ly: " << Ly << endl;
    fout << "seed: " << seed << endl;
    fout << "np: " << tot_rank << endl;
  }
}

void BasicDomain::create_cell_list() {
  for (int i = 0; i < end_pos; i++) {
    int idx = particle[i].cell_idx(yl, ncols);
    cell[idx].push_front(&particle[i]);
  }
}

void BasicDomain::ini_disorder(double *disorder, int n) {
  double d = 1.0 / (n - 1);
  for (int i = 0; i < n; i++)
    disorder[i] = (-0.5 + i * d) * eps * 2.0 * PI;
  shuffle(disorder, n, myran);
}

void BasicDomain::update_position_edge_row(int row) {
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
        if (!disorder_free) {
          noise += cell[curNode->cell_idx(yl, ncols)].disorder;
        }
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

void BasicDomain::update_position_inner_rows() {
  for (int row = 1; row < nrows - 1; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      cell[col + j].head = NULL;
      cell[col + j].size = 0;
    }
  }
  my_nPar = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!particle[i].is_empty) {
      if (!particle[i].is_moved) {
        double noise = eta * 2.0 * PI * (myran->doub() - 0.5);
        if (!disorder_free) {
          noise += cell[particle[i].cell_idx(yl, ncols)].disorder;
        }
        particle[i].update_coor(noise, Lx, yl);
        particle[i].is_ghost = false;
      }
      if (!particle[i].is_ghost) {
        int idx = particle[i].cell_idx(yl, ncols);
        cell[idx].push_front(&particle[i]);
        my_nPar++;
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
  int nPar = buf_size / 4;
  bool ghost;
  if (row == 0 || row == nrows - 1) {
    ghost = true;
  } else if (row > 0 && row < nrows - 1) {
    ghost = false;
    my_nPar += nPar;
  } else {
    cout << "Failed to unpack, row = " << row << endl;
    exit(1);
  }
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

void BasicDomain::update_velocity() {
  MPI_Request req[4];
  MPI_Status stat[4];
  double *buf[4];
  int buf_size[4];
  for (int i = 0; i < 4; i++) {
    buf[i] = new double[MAX_BUF_SIZE];
    buf_size[i] = MAX_BUF_SIZE;
  }

  /* pre_rank -> myrank */
  recv(buf[0], buf_size[0], pre_rank, 13, &req[0]);
  /* myrank -> next_rank */
  send(nrows - 2, buf[1], buf_size[1], next_rank, 13, &req[1]);
  /* myrank <- next_rank */
  recv(buf[2], buf_size[2], next_rank, 31, &req[2]);
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
void BasicDomain::shift_boundary(int pre_nPar, int next_nPar, int * offset) {
  int threshold = 2 * nPar_per_row;
  if (myrank == 0) {
    offset[0] = 0;
    int dn = next_nPar - my_nPar;
    if (dn > threshold) {
      offset[1] = 1;
    } else if (dn < -threshold) {
      offset[1] = -1;
    } else {
      offset[1] = 0;
    }
  } else if (myrank == tot_rank - 1) {
    offset[1] = 0;
    int dn = pre_nPar - my_nPar;
    if (dn > threshold) {
      offset[0] = 1;
    } else if (dn < -threshold) {
      offset[0] = -1;
    } else {
      offset[0] = 0;
    }
  } else {
    int dn = pre_nPar - my_nPar;
    if (dn > threshold) {
      offset[0] = 1;
    } else if (dn < -threshold) {
      offset[0] = -1;
    } else {
      offset[0] = 0;
    }
    dn = next_nPar - my_nPar;
    if (dn > threshold) {
      offset[1] = 1;
    } else if (dn < -threshold) {
      offset[1] = -1;
    } else {
      offset[1] = 0;
    }
  }
}
void BasicDomain::update_position() {
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
  update_position_edge_row(nrows - 2);
  send(nrows - 1, buf[1], buf_size[1], next_rank, 12, &req[1]);
  /* myrank <- next_rank */
  recv(buf[2], buf_size[2], next_rank, 21, &req[2]);
  /* pre_rank <- myrank */
  update_position_edge_row(1);
  send(0, buf[3], buf_size[3], pre_rank, 21, &req[3]);

  /* update position of inner rows */
  update_position_inner_rows();

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
StaticDomain::StaticDomain(const cmdline::parser & cmd): BasicDomain(cmd) {
  set_cell_list(eps);
  ini_output("s", eta, eps, rho0, cmd.get<unsigned long long>("seed"));
}

StaticDomain::~StaticDomain() {
  delete[] cell;
  delete[] particle;
  if (myrank == 0) {
    fout.close();
  } else if (myrank == 1) {
    fout << "-------- End --------\n";
    fout << "Finished at " << asctime(localtime(&beg_time));
    fout.close();
  }
  if (cgs)
    delete cgs;
}

void StaticDomain::set_cell_list(double eps) {
  cell = new Cell[ncols * nrows];
  Cell::find_all_neighbor(cell, ncols, nrows);
  if (eps > 0) {
    int my_n = ncols * (nrows - 2);
    int n = my_n * tot_rank;
    double *my_disorder = new double[my_n];
    double *disorder = new double[n];
    if (myrank == 0) {
      ini_disorder(disorder, n);
    }
    MPI_Scatter(disorder, my_n, MPI_DOUBLE,
                my_disorder, my_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int row = 1; row < nrows - 1; row++) {
      for (int col = 0; col < ncols; col++) {
        int i = col + row * ncols;
        int j = col + (row - 1) * ncols;
        cell[i].disorder = my_disorder[j];
      }
    }
    delete[] disorder;
    delete[] my_disorder;
  }
  double sum = 0;
  for (int row = 1; row < nrows - 1; row++) {
    for (int col = 0; col < ncols; col++) {
      int i = col + row * ncols;
      sum += cell[i].disorder;
    }
  }
  double tot_sum = 0;
  MPI_Reduce(&sum, &tot_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    cout << "sum disorder = " << tot_sum << endl;
  }
}

void StaticDomain::one_step(int t) {
  update_velocity();
  update_position();
  output(t);
}

/****************************************************************************/
/*                   Dynamic domain decomposition                           */
/****************************************************************************/
DynamicDomain::DynamicDomain(const cmdline::parser & cmd): BasicDomain(cmd) {
  set_cell_list();
  row_offset[0] = row_offset[1] = 0;
  refresh_rate = cmd.get<int>("rate");
  char buf[10];
  snprintf(buf, 10, "d%d", refresh_rate);
  ini_output(buf, eta, eps, rho0, cmd.get<unsigned long long>("seed"));
}

DynamicDomain::~DynamicDomain() {
  delete[] cell_buf;
  delete[] particle;
  if (myrank == 0) {
    fout.close();
  } else if (myrank == 1) {
    fout << "-------- End --------\n";
    fout << "Finished at " << asctime(localtime(&beg_time));
    fout.close();
  }
  if (cgs)
    delete cgs;
}

void DynamicDomain::set_cell_list() {
  int row_excess = (nrows - 2) / 2;
  int row_beg, row_end;
  if (myrank == 0) {
    MAX_ROW = nrows + row_excess;
    CELL_BUF_SIZE = ncols * MAX_ROW;
    cell_buf = new Cell[CELL_BUF_SIZE];
    cell = cell_buf;
    row_beg = -1;
    row_end = MAX_ROW - 1;
  } else if (myrank == tot_rank - 1) {
    MAX_ROW = nrows + row_excess;
    CELL_BUF_SIZE = ncols * MAX_ROW;
    cell_buf = new Cell[CELL_BUF_SIZE];
    cell = cell_buf + row_excess * ncols;

    row_end = (nrows - 2) * tot_rank;
  } else {
    MAX_ROW = nrows + 2 * row_excess;
    CELL_BUF_SIZE = ncols * MAX_ROW;
    cell_buf = new Cell[CELL_BUF_SIZE];
    cell = cell_buf + row_excess * ncols;
  }
  Cell::find_all_neighbor(cell, ncols, nrows);

  if (eps > 0) {
    int tot_cell = ncols * (nrows - 2) * tot_rank;
    double *disorder = new double[tot_cell];
    if (myrank == 0) {
      ini_disorder(disorder, tot_cell);
      double sum = 0;
      for (int i = 0; i < tot_cell; i++)
        sum += disorder[i];
      cout << "sum of disorder = " << sum << endl;
    }
    MPI_Bcast(disorder, tot_cell, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      for (int row = 1; row < MAX_ROW; row++) {
        for (int col = 0; col < ncols; col++) {
          int i = col + row * ncols;
          int j = col + (row - 1) * ncols;
          cell_buf[i].disorder = disorder[j];
        }
      }
    } else if (myrank < tot_rank - 1){
      for (int row = 0; row < MAX_ROW; row++) {
        for (int col = 0; col < ncols; col++) {
          int i = col + row * ncols;
          int j = col + ((nrows - 2) * myrank - 1 - row_excess + row ) * ncols;
          cell_buf[i].disorder = disorder[j];
        }
      }
    } else {
      for (int row = 0; row < MAX_ROW - 1; row++) {
        for (int col = 0; col < ncols; col++) {
          int i = col + row * ncols;
          int j = col + ((nrows - 2) * myrank - 1 - row_excess + row) * ncols;
          cell_buf[i].disorder = disorder[j];
        }
      }
    }
    delete[] disorder;
  }
  {
    double sum = 0;
    for (int row = 1; row < nrows - 1; row++) {
      for (int col = 0; col < ncols; col++) {
        int i = col + row * ncols;
        sum += cell[i].disorder;
      }
    }
    double tot_sum;
    MPI_Reduce(&sum, &tot_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0)
      cout << "sum of disorder = " << tot_sum << endl;
  }
}

void DynamicDomain::global_rearrange(int t) {
  if (t % refresh_rate == 0) {
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
  } else {
    row_offset[0] = row_offset[1] = 0;
  }
}

void DynamicDomain::local_rearrange(int t) {
  if (t % refresh_rate == 0) {
    int pre_nPar, next_nPar;
    int pre_id, next_id;
    if (myrank == 0) {
      pre_id = MPI_PROC_NULL;
      next_id = next_rank;
    } else if (myrank == tot_rank - 1) {
      pre_id = pre_rank;
      next_id = MPI_PROC_NULL;
    } else {
      pre_id = pre_rank;
      next_id = next_rank;
    }
    MPI_Sendrecv(&my_nPar, 1, MPI_INT, next_id, 45,
                 &pre_nPar, 1, MPI_INT, pre_id, 45,
                 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    MPI_Sendrecv(&my_nPar, 1, MPI_INT, pre_id, 54,
                 &next_nPar, 1, MPI_INT, next_id, 54,
                 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    shift_boundary(pre_nPar, next_nPar, row_offset); 
  }
}

void DynamicDomain::update_velocity_dynamic() {
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

  MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
  for (int i = 0; i < 4; i++) {
    delete[] buf[i];
    buf[i] = NULL;
  }
}

void DynamicDomain::one_step(int t) {
  //global_rearrange(t);
  local_rearrange(t);
  update_velocity_dynamic();
  update_position();
  output(t);
}
