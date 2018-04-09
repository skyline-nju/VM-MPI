#ifndef COMN_H
#define COMN_H
//#define USE_MPI

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>
#include <ctime>
#include <chrono>
#include <algorithm>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

const double PI = 3.14159265358979;

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

// create folder
void mkdir(const char *folder);

// split string by a delimiter
std::vector<std::string> split(
  const std::string &str, const std::string &dlm);

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

// Calculate packing fraction from particle number
double cal_packing_fraction_2(int n, double Lx, double Ly, double sigma);

// Calculate particle number from packing fraction
int cal_particle_number_2(double phi, double Lx, double Ly, double sigma);

// Calculate the elapsed time to run a function
template <typename _TFunc>
void cal_elapsed_time(_TFunc f) {
  auto t1 = std::chrono::system_clock::now();
  f();
  auto t2 = std::chrono::system_clock::now();
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

