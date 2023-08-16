#ifndef COMN_H
#define COMN_H
//#define USE_MPI

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include "vect.h"

/*************************************************************************//**
 *                      Global constans
 *****************************************************************************/

const double PI = 3.14159265358979323846; //!< pi

//!< path delimiter, "\\" for windows and '/' for linux
#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

/*************************************************************************//**
 *
 *                      Periodic boundary condition
 *                      
 *****************************************************************************/

/*
 * \brief Shift the x if out of the range.
 * \tparam T1     Template type 1
 * \tparam T2     Template type 2
 * \param x       Input/output coordinate
 * \param x_min   Min of x
 * \param x_max   Max of x
 * \param len     Distance from x_min to x_max
 */
template <typename T1, typename T2>
void tangle_1(T1 &x, T2 x_min, T2 x_max, T2 len) {
  if (x < x_min) {
    x += len;
  } else if (x >= x_max) {
    x -= len;
  }
}

template <typename T1, typename T2>
void tangle_1(T1 &x, T2 len) {
  if (x < 0) {
    x += len;
  } else if (x >= len) {
    x -= len;
  }
}

template <typename T1, typename T2>
void tangle_3(Vec_3<T1> &x, const Vec_3<T2> &x_min, const Vec_3<T2> &x_max, const Vec_3<T2> &len) {
  tangle_1(x.x, x_min.x, x_max.x, len.x);
  tangle_1(x.y, x_min.y, x_max.y, len.y);
  tangle_1(x.z, x_min.z, x_max.z, len.z);
}

template <typename T1, typename T2>
void tangle_3(Vec_3<T1> &x, const Vec_3<T2> &len) {
  tangle_1(x.x, len.x);
  tangle_1(x.y, len.y);
  tangle_1(x.z, len.z);
}

template <typename T1, typename T2>
void tangle_3(Vec_3<T1> &x, const Vec_3<T2> &x_min, const Vec_3<T2> &x_max, const Vec_3<T2> &len, const Vec_3<bool> &flag_ext) {
  if (!flag_ext.x) tangle_1(x.x, x_min.x, x_max.x, len.x);
  if (!flag_ext.y) tangle_1(x.y, x_min.y, x_max.y, len.y);
  if (!flag_ext.z) tangle_1(x.z, x_min.z, x_max.z, len.z);
}

template <typename T1, typename T2>
void tangle_3(Vec_3<T1> &x, const Vec_3<T2> &len, const Vec_3<bool> &flag_ext) {
  if (!flag_ext.x) tangle_1(x.x, len.x);
  if (!flag_ext.y) tangle_1(x.y, len.y);
  if (!flag_ext.z) tangle_1(x.z, len.z);
}
/**
 * \brief Cal nearest distance under the periodic boundary condition
 * \tparam T       Template type
 * \param dx       Input/output distance
 * \param len      Length of the domain 
 * \param half_len Half length of the domain
 */
template <typename T>
void untangle_1(T &dx, T len, T half_len) {
  if (dx < -half_len) {
    dx += len;
  } else if (dx > half_len) {
    dx -= len;
  }
}

template <typename T>
void untangle_3(Vec_3<T> &dis, const Vec_3<T> &len, const Vec_3<T> &half_len) {
  untangle_1(dis.x, len.x, half_len.x);
  untangle_1(dis.y, len.y, half_len.y);
  untangle_1(dis.z, len.z, half_len.z);
}

template <typename T>
void untangle_3(Vec_3<T> &dis, const Vec_3<T> &len, const Vec_3<T> &half_len, const Vec_3<bool> &flag_ext) {
  if (!flag_ext.x) untangle_1(dis.x, len.x, half_len.x);
  if (!flag_ext.y) untangle_1(dis.y, len.y, half_len.y);
  if (!flag_ext.z) untangle_1(dis.z, len.z, half_len.z);
}

/*************************************************************************//**
 * \brief Base class for data exporter.
 ****************************************************************************/
class BaseExporter {
public:
  BaseExporter(): n_step_(0), frame_interval_(0), iframe_(0) {}
  BaseExporter(int n_step, int sep, int start = 0);
  void set_lin_frame(int sep, int n_step, int start  = 0);
  bool need_export(int i_step);
  size_t get_n_frames() const { return frames_arr_.size(); }
protected:
  int n_step_;
  int frame_interval_;
private:
  int iframe_;
  std::vector<int> frames_arr_;
};

/*************************************************************************//**
 * \brief Base class for log exporter
 * 
 * To output parameters, the constructor of the derived class should contains
 * some lines:
 * fout_ << "\n-------- Parameters --------";
 * fout_ << "\npara_name=" << para_value;
 * fout_ << "\n\n-------- RUN --------";
 * fout_ << "\ntime step\telapsed time" << std::endl;
 ****************************************************************************/
class BaseLogExporter:public BaseExporter {
public:
  explicit BaseLogExporter(const std::string &filename, int n_par,
                           int end, int sep, int start=0);

  ~BaseLogExporter();

  void record(int i_step);
protected:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  std::ofstream fout_;
private:
  int n_par_;
};

/*************************************************************************//**
 * \brief Create a folder
 * \param folder The name of the folder to be created
 ****************************************************************************/
void mkdir(const char *folder);

inline void mkdir(const std::string &folder) {mkdir(folder.c_str());}

// split string by a delimiter
std::vector<std::string> split(const std::string &str, const std::string &dlm);

template <class T>
void str_to_num(const std::string str, T &num) {
  std::stringstream ss;
  ss << str;
  ss >> num;
}

template <class T>
void num_to_str(const T &num, std::string str) {
  std::stringstream ss;
  ss << num;
  ss >> str;
}

// Calculate packing fraction from particle number in 2D domain
double cal_packing_fraction_2(int n, double Lx, double Ly, double sigma);

// Calculate particle number from packing fraction in 2D domain
int cal_particle_number_2(double phi, double Lx, double Ly, double sigma);

// Calculate the elapsed time to run a function
template <typename TFunc>
void cal_elapsed_time(TFunc f) {
  const auto t1 = std::chrono::system_clock::now();
  f();
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}

//#define SPATIAL_SORT
#ifdef SPATIAL_SORT
#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort.h>
struct MyLessX {
  template <class MyPoint>
  bool operator()(const MyPoint &p, const MyPoint &q) const {
    return p.x < q.x;
  }
};

struct MyLessY {
  template <class MyPoint>
  bool operator()(const MyPoint &p, const MyPoint &q) const {
    return p.y < q.y;
  }
};

template <class MyPoint>
struct MySpatialSortingTraits {
  typedef MyPoint Point_2;
  typedef MyLessX Less_x_2;
  typedef MyLessY Less_y_2;

  Less_x_2 less_x_2_object() const {
    return Less_x_2();
  }
  Less_y_2 less_y_2_object() const {
    return Less_y_2();
  }
};
#endif

#endif

