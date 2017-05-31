#include "node.h"
#include "mpi.h"

using namespace std;
double Node::v0 = 0.5;

void Node::align(Node *node) {
  if (rr(node) < 1)
    addV(node);
}

void Node::align(Node *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

void Node::update_coor(double noise, double Lx, double yl) {
  double tmp = sqrt(vx*vx + vy*vy);
  double c1 = vx / tmp;
  double s1 = vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx0 = c1 * c2 - s1 * s2;
  vy = vy0 = c1 * s2 + c2 * s1;
  x += v0*vx;
  if (x >= Lx) {
    x -= Lx;
  } else if (x < 0) {
    x += Lx;
  }
  y += v0*vy;
  is_moved = true;
  new_arrival = false;
}

void Node::ini_random(Node * par, int nPar, Ran * myran,
                      double Lx, double yl, double yh) {
  double Ly0 = yl + 1;
  double Ly1 = yh - 1;
  for (int i = 0; i < nPar; i++) {
    par[i].x = myran->doub() * Lx;
    par[i].y = myran->doub() * (Ly1 - Ly0) + Ly0;
    double theta = myran->doub() * 2 * PI;
    par[i].vx = par[i].vx0 = cos(theta);
    par[i].vy = par[i].vy0 = sin(theta);
    par[i].is_moved = false;
    par[i].is_empty = false;
    par[i].is_ghost = false;
    par[i].new_arrival = false;
  }
}

void Node::ini_from_snap(Node **par, int &npar, int &max_par_num,
                         const std::string filename, double Lx, double Ly,
                         int tot_rank, int myrank) {
  MPI_File fin;
  if (!exist(filename.c_str())) {
    cout << "Error, " << filename << " doesn't exist!\n";
    exit(1);
  }
  MPI_File_open(
    MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
  MPI_Offset size;
  MPI_File_get_size(fin, &size);
  int size_buff = size / 4;
  float *buff = new float[size_buff];
  MPI_Status stat;
  MPI_File_read_at(fin, 0, buff, size_buff, MPI_FLOAT, &stat);

  vector<string> str_list = split(filename, "_");
  double lx, ly;
  str_to_num(str_list[4], lx);
  str_to_num(str_list[5], ly);
  int nx, ny;
  if (int(Lx) % int(lx) == 0 && int(Ly) / tot_rank % int(ly) == 0) {
    nx = int(Lx) / int(lx);
    ny = int(Ly) / tot_rank / int(ly);
  } else {
    cout << "Error, the size of input snapshot is not right" << endl;
    exit(1);
  }
  int npar_in = size_buff / 3;
  int idx = 0;
  max_par_num = npar * nx * ny * 3;
  Node *particle = new Node[max_par_num];
  for (int row = 0; row < ny; row++) {
    double dy = row * ly + Ly / tot_rank * myrank;
    for (int col = 0; col < nx; col++) {
      double dx = col * lx;
      for (int i = 0; i < npar_in; i++) {
        particle[idx].x = buff[i * 3] + dx;
        particle[idx].y = buff[i * 3 + 1] + dy;
        particle[idx].vx = particle[idx].vx0 = cos(buff[i * 3 + 2]);
        particle[idx].vy = particle[idx].vy0 = sin(buff[i * 3 + 2]);
        particle[i].is_moved = false;
        particle[i].is_empty = false;
        particle[i].is_ghost = false;
        particle[i].new_arrival = false;
        idx++;
      }
    }
  }
  npar = idx;
  *par = particle;
  delete[] buff;
}

