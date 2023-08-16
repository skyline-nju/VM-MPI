#include "comn.h"
#include <ctime>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <chrono>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif


using namespace std;

void mkdir(const char *folder) {
#ifdef _MSC_VER
  if (_access(folder, 0) != 0)
#else
  if (access(folder, 0) != 0)
#endif
  {
    char command[100];
    snprintf(command, 100, "mkdir -p %s", folder);
    if (system(command))
      cout << "create folder: " << folder << " successfully" << endl;
  } else
    cout << "folder: " << folder << " already exists" << endl;
}

vector<string> split(const string &str, const string &delim) {
  string::size_type pos;
  vector<string> res;
  int str_size = str.size();
  int dlm_size = delim.size();
  for (unsigned int i = 0; i < str.size(); i++) {
    pos = str.find(delim, i);
    if (pos == string::npos) {
      res.push_back(str.substr(i, str_size));
      break;
    } else {
      res.push_back(str.substr(i, pos - i));
      i = pos + dlm_size - 1;
    }
  }
  return res;
}

double cal_packing_fraction_2(int n, double Lx, double Ly, double sigma) {
  double a = sigma * 0.5;
  return PI * a * a * n / (Lx * Ly);
}

int cal_particle_number_2(double phi, double Lx, double Ly, double sigma) {
  double phi_max = PI / (2 * sqrt(3.0));
  if (phi > phi_max) {
    cout << "Input packing fraction phi = " << phi
      << " is larger than phi_max = " << phi_max << endl;
    exit(1);
  } else {
    double a = sigma * 0.5;
    return int(round(phi * Lx * Ly / (PI * a * a)));
  }
}
/**************************************************************************//**
 * @brief Construct a new Base Exporter:: Base Exporter object
 * 
 * @param n_step Total time steps to run
 * @param sep    Time seperation to dump frame.
 * @param start  Starting time step to dump frame.
 ****************************************************************************/
BaseExporter::BaseExporter(int n_step, int sep, int start)
  : n_step_(n_step), frame_interval_(sep), iframe_(0) {
  frames_arr_.reserve((n_step - start) / sep);
  for (auto i = start + sep; i <= n_step_; i += sep) {
    frames_arr_.push_back(i);
  }
}

void BaseExporter::set_lin_frame(int sep, int n_step, int start) {
  if (sep <= 0) {
    std::cerr << "line " << __LINE__ << " of " << __FILE__
      << ": sep must be positive!" << std::endl;
    exit(1);
  }
  frame_interval_ = sep;
  n_step_ = n_step;
  for (auto i = start + sep; i <= n_step_; i += sep) {
    frames_arr_.push_back(i);
  }
}

bool BaseExporter::need_export(int i_step) {
  auto flag = false;
  if (!frames_arr_.empty() && i_step == frames_arr_[iframe_]) {
    iframe_++;
    flag = true;
  }
  return flag;
}
/*************************************************************************//**
 * @brief Construct a new Base Log Exporter:: Base Log Exporter object
 * 
 * @param filename    Filename of log file.
 * @param n_par       Num of particles
 * @param n_step      Total time steps to run.
 * @param sep         Frame spacing
 * @param start       First time step to dump frame.
 ***************************************************************************/
BaseLogExporter::BaseLogExporter(const std::string& filename, int n_par,
                                 int n_step, int sep, int start)
  : BaseExporter(n_step, sep, start), fout_(filename), n_par_(n_par) {
  t_start_ = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(t_start_);
  char str[100];
  // ReSharper disable CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&start_time));
  // ReSharper restore CppDeprecatedEntity
  fout_ << "Started simulation at " << str << "\n";
}

BaseLogExporter::~BaseLogExporter() {
  const auto t_now = std::chrono::system_clock::now();
  auto end_time = std::chrono::system_clock::to_time_t(t_now);
  char str[100];
  // ReSharper disable CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&end_time));
  // ReSharper restore CppDeprecatedEntity
  fout_ << "Finished simulation at " << str << "\n";
  std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
  fout_ << "speed=" << std::scientific << n_step_ * double(n_par_) / elapsed_seconds.count()
    << " particle time step per second per core\n";
  fout_.close();
}

void BaseLogExporter::record(int i_step) {
  if (need_export(i_step)) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = int(dt / 3600);
    const auto min = int((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout_ << i_step << "\t" << hour << ":" << min << ":" << sec << std::endl;
  } 
}
