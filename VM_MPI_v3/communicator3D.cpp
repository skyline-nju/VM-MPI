#include "communicator3D.h"
#include "comn.h"
#include "rand.h"

/*
void test_comm_velocity() {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int max = 1000000;
  double *data = new double[max]{};
  Ran myran(1);
  if (my_rank == 0) {
    for (int i = 0; i < max; i++) {
      data[i] = myran.doub();
    }  
  }
  int n = 1000;
  int *idx_arr = new int[n]{};
  if (my_rank == 0) {
    for (int i = 0; i < n; i++) {
      idx_arr[i] = int(myran.doub() * max);
    }
  }

  auto comm1 = [n, &idx_arr, &data, my_rank]() {
    for (int i = 0; i < n; i++) {
      if (my_rank == 0) {
        MPI_Send(&data[idx_arr[i]], 1, MPI_DOUBLE, 1, i + 1, MPI_COMM_WORLD);
      } else if (my_rank == 1) {
        MPI_Recv(&data[i], 1, MPI_DOUBLE, 0, i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  };

  auto comm2 = [n, &idx_arr, &data, my_rank]() {
    double *buff = new double[n];
    if (my_rank == 0) {
      for (int i = 0; i < n; i++) {
        buff[i] = data[idx_arr[i]];
      }
    }
    if (my_rank == 0) {
      MPI_Send(buff, n, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
    } else if (my_rank == 1) {
      MPI_Recv(buff, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (my_rank == 1) {
      for (int i = 0; i <n; i++) {
        data[i] = buff[i];
      }
    }
    delete[] buff;
  };

  auto commI2 = [n, &idx_arr, &data, my_rank]() {
    double *buff = new double[n];
    if (my_rank == 0) {
      for (int i = 0; i < n; i++) {
        buff[i] = data[idx_arr[i]];
      }
    }

    MPI_Request req;
    MPI_Status stat;
    if (my_rank == 0) {
      MPI_Isend(buff, n, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &req);
    } else if (my_rank == 1) {
      MPI_Irecv(buff, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req);
    }

    MPI_Wait(&req, &stat);
    if (my_rank == 1) {
      for (int i = 0; i < n; i++) {
        data[i] = buff[i];
      }
    }
    delete[] buff;
  };

  auto comm3 = [n, &idx_arr, &data, my_rank]() {
    int pos;
    int buff_size = sizeof(double) * n + 100;
    char *buff = new char[buff_size];
    if (my_rank == 0) {
      for (int i = 0; i < n; i++) {
        MPI_Pack(&data[idx_arr[i]], 1, MPI_DOUBLE, buff, buff_size, &pos, MPI_COMM_WORLD);
      }
    }
    if (my_rank == 0) {
      MPI_Send(buff, buff_size, MPI_PACKED, 1, 1, MPI_COMM_WORLD);
    } else if (my_rank == 1) {
      MPI_Recv(buff, buff_size, MPI_PACKED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (my_rank == 1) {
      for (int i = 0; i < n; i++) {
        MPI_Unpack(buff, buff_size, &pos, &data[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
    delete[] buff;
  };

  auto show = [my_rank, &data, &idx_arr](int i) {
    if (my_rank == 0) {
      std::cout << data[idx_arr[i]] << std::endl;
    } else if (my_rank == 1) {
      std::cout << data[i] << std::endl;
    }
  };

  //MPI_Barrier(MPI_COMM_WORLD);
  //cal_elapsed_time(comm1);
  //show(199);

  MPI_Barrier(MPI_COMM_WORLD);
  cal_elapsed_time(comm2);
  show(199);

  //MPI_Barrier(MPI_COMM_WORLD);
  //cal_elapsed_time(commI2);
  //show(199);

  //MPI_Barrier(MPI_COMM_WORLD);
  //cal_elapsed_time(comm3);
  //show(199);

}
*/

void find_shell(const Vec_3<int> &n, const Vec_3<int> &thickness,
                Vec_3<block_t> shell[2]) {
  for (int ori = 0; ori < 3; ori++) {
    if (thickness[ori]) {
      for (int dim = 0; dim < 3; dim++) {
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

void set_comm_block(const Vec_3<int> &cells_size, const Vec_3<bool> &flag_comm,
                    Vec_3<block_t> inner_shell[2], Vec_3<block_t> outer_shell[2],
                    MPI_Comm group_comm) {
  Vec_3<int> thickness{};
  for (int dim = 0; dim < 3; dim++) {
    thickness[dim] = flag_comm[dim] ? 1 : 0;
  }

  find_shell(cells_size, -thickness, inner_shell);
  for (int ori = 0; ori < 3; ori++) {
    inner_shell[0][ori].beg += thickness;
    inner_shell[0][ori].end += thickness;
    inner_shell[1][ori].beg += thickness;
    inner_shell[1][ori].end += thickness;
  }

  const Vec_3<int> extended_cells_size = cells_size + thickness * 2;
  find_shell(extended_cells_size, thickness, outer_shell);


  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  if (my_rank == 0) {
    std::cout << "proc = " << my_rank << std::endl;
    std::cout << "n: " << cells_size << std::endl;
    std::cout << "thickness: " << thickness << std::endl;

    std::cout << "outer shell" << std::endl;
    for (int i = 0; i < 2; i++) {
      for (int dim = 0; dim < 3; dim++) {
        std::cout << outer_shell[i][dim].beg << ";\t" << outer_shell[i][dim].end << "\n";
      }
    }  
    std::cout << "inner shell" << std::endl;
    for (int i = 0; i < 2; i++) {
      for (int dim = 0; dim < 3; dim++)
        std::cout << inner_shell[i][dim].beg << ";\t" << inner_shell[i][dim].end << "\n";  
    }
    std::cout << std::endl;
  std::cout << "----------" << std::endl;
  }
}
