#include "subDomain.h"

DomainBase::DomainBase(double Lx, double Ly, double Rho, double eta0,
                       Ull seed0, double v)
    : system_size(Lx, Ly), rho0(Rho), eta(eta0), seed(seed0), v0(v) {
  nPar = int(Lx * Ly * Rho);
  std::cout << "nPar = " << nPar << "\n";
}

DomainBase::~DomainBase() {
  delete myran;
}

#ifndef USE_MPI
SingleDomain::SingleDomain(double Lx, double Ly, double Rho, double eta0, Ull seed0)
  : DomainBase(Lx, Ly, Rho, eta0, seed0) {
  pbc = new PBC_2(Lx, Ly);
  myran = new Ran(seed);
}

void SingleDomain::ini_rand() {
  p_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    p_arr.emplace_back(*myran, system_size.x, system_size.y);
  }
  cl = new CellListIdx_2(system_size.x, system_size.y, 1);
  cl->create(p_arr);
}

void SingleDomain::cal_force() {
  auto f_ij = [this](int i, int j) {
    p_arr[i].interact(p_arr[j], *pbc);
  };
  cl->for_each_pair(f_ij, 0);
}

void SingleDomain::integrate() {
  for (int i = 0; i < nPar; i++) {
    p_arr[i].move(eta, v0, *myran, *pbc);
  }
  cl->update(p_arr);
}

void SingleDomain::cal_order_para(double & phi, double & theta) const {
  double svx = 0;
  double svy = 0;
  for (int i = 0; i < nPar; i++) {
    svx += p_arr[i].vx;
    svy += p_arr[i].vy;
  }
  phi = std::sqrt(svx * svx + svy * svy) / nPar;
  theta = std::atan2(svy, svx);
}


SingleDomain2::SingleDomain2(double Lx, double Ly, double Rho, double eta0, Ull seed0)
  : DomainBase(Lx, Ly, Rho, eta0, seed0) {
  pbc = new PBC_2(Lx, Ly);
  myran = new Ran(seed);

  std::cout << "seed = " << seed << "\n";
}

void SingleDomain2::ini_rand() {
  p_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    p_arr.emplace_back(*myran, system_size.x, system_size.y);
  }
  cl = new CellListNode_2(system_size.x, system_size.y, 1.0);
  cl->create(p_arr);

  double phi, theta;
  cal_order_para(phi, theta);
  std::cout << "t = 0\t phi = " << phi << ", theta = " << theta << "\n";

}

void SingleDomain2::cal_force() {
  auto f_ij = [this](Node *pi, Node *pj) {
    pi->interact(*pj, *pbc);
  };
  cl->for_each_pair(f_ij, 0);
}

void SingleDomain2::integrate() {
  for (int i = 0; i < nPar; i++) {
    p_arr[i].move(eta, v0, *myran, *pbc);
  }
  cl->recreate(p_arr);
}

void SingleDomain2::integrate2() {
  auto lambda = [this](Node *p) {
    p->move(eta, v0, *myran, *pbc);
  };
  cl->for_each(lambda);
  cl->recreate(p_arr);
}

void SingleDomain2::cal_order_para(double & phi, double & theta) const {
  func_order_para(p_arr, phi, theta);
}

SingleDomain3::SingleDomain3(double Lx, double Ly, double Rho, double eta0, Ull seed0)
  : DomainBase(Lx, Ly, Rho, eta0, seed0) {
  pbc = new PBC_2(Lx, Ly);
  myran = new Ran(seed);
  std::cout << "v0 = " << v0 << std::endl;
}

void SingleDomain3::ini_rand() {
  p_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    p_arr.emplace_back(*myran, system_size.x, system_size.y);
  }
  cl = new CellListBiNode_2(system_size.x, system_size.y, 1.0);
  cl->create(p_arr);
}

void SingleDomain3::cal_force() {
  auto f_ij = [this](BiNode *pi, BiNode *pj) {
    pi->interact(*pj, *pbc);
  };
  cl->for_each_pair(f_ij, 0);
}

void SingleDomain3::integrate() {
  //for (int i = 0; i < nPar; i++) {
  //  p_arr[i].move(eta, v0, *myran, *pbc);
  //}
  //cl->recreate(p_arr);
  auto lambda = [this](BiNode *p) {
    p->move(eta, v0, *myran, *pbc);
  };
  cl->update(lambda);
}

void SingleDomain3::integrate2() {

  auto lambda = [this](BiNode *p, Vec_2<int> &d_cell) {
    p->move(eta, v0, *myran, *pbc, d_cell);
  };
  cl->update_by_row(lambda);
  //cl->update2(lambda);
}

void SingleDomain3::cal_order_para(double & phi, double & theta) const {
  double svx = 0;
  double svy = 0;
  for (int i = 0; i < nPar; i++) {
    svx += p_arr[i].vx;
    svy += p_arr[i].vy;
  }
  phi = std::sqrt(svx * svx + svy * svy) / nPar;
  theta = std::atan2(svy, svx);
}


#else
StaticSubDomain::StaticSubDomain(double Lx, double Ly, double Rho, double eta0, Ull seed0)
  : DomainBase(Lx, Ly, Rho, eta0, seed0) {
  MPI_Comm_size(MPI_COMM_WORLD, &tot_processor);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  pre_rank = my_rank == 0 ? tot_processor - 1 : my_rank - 1;
  next_rank = my_rank == tot_processor - 1 ? 0 : my_rank + 1;
  double Ly_sub = Ly / tot_processor;
  cl = new CellListIdx_2(Lx, Ly_sub, Ly_sub * my_rank, 1.0);
  pbc_x = new PBC_X_2(Lx);
  myran = new Ran(seed + my_rank);
  my_nPar = nPar / tot_processor;
}

void StaticSubDomain::ini_rand() {
  p_arr.reserve(1.5 * my_nPar);
  double y0 = system_size.y / tot_processor * my_rank;
  for (int i = 0; i < my_nPar; i++) {
    p_arr.emplace_back(*myran, system_size.x, system_size.y, 0., y0);
  }
  cl->create(p_arr);
}

void StaticSubDomain::cal_force() {
  auto f_ij = [this](int i, int j) {
    p_arr[i].interact(p_arr[j], *pbc_x);
  };


}

void StaticSubDomain::integrate() {
}

#endif


