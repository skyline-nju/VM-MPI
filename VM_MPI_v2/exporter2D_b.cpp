#include "exporter2D_b.h"

SnapExporter::SnapExporter(int frame_interval, int n_step): 
  BaseExporter(n_step, frame_interval) {
  std::string path = "data" + delimiter + "snap" + delimiter;
  create_output_folder(path);
  snprintf(file_prefix_, 100, "%ssnap_%s", path.c_str(), set_base_name().c_str());
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc_);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
}

RhoxExpoerter::RhoxExpoerter(int sep, int n_step, const Grid_2& grid, const Vec_2<double>& origin):
  BaseExporter(n_step, sep), origin_(origin) {
  std::string path = "data" + delimiter;
  create_output_folder(path);
  char filename[100];
  snprintf(filename, 100, "%srhoxUINT16_%s.bin", path.c_str(), set_base_name().c_str());
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc_);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

  frame_size_ = grid.gl_n().x * sizeof(unsigned short);
  offset_ = grid.origin().x * sizeof(unsigned short);

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, 
                MPI_INFO_NULL, &fh_);
  buf_size_ = grid.n().x;
  buf_ = new unsigned short[buf_size_];
}
