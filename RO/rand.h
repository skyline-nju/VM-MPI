#ifndef RAND_H
#define RAND_H
#include <cmath>
#include <vector>
/********************Uniform Random number generator*************************/
// Ref: Numerical Recipes, The Art of Scientific Computing, 3rd
// Uniform Random number generator with period 3.138e57
struct Ran {
  unsigned long long u, v, w;
  Ran(unsigned long long j) :v(4101842887655102017LL), w(1) {
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }
  inline unsigned long long int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }
  // Returns random double-precision floating value between 0. and 1.
  inline double doub() {
    return 5.42101086242752217E-20 * int64();
  }
  inline unsigned int int32() {
    return (unsigned int)int64();
  }
};


// Fastest Uniform Random number generator with period 1.8e19
struct Ranq1 {
  unsigned long long v;
  Ranq1(unsigned long long j) : v(4101842887655102017LL) {
    v ^= j;
    v = int64();
  }
  inline unsigned long long int64() {
    v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
    return v * 2685821657736338717LL;
  }
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  inline unsigned int int32() { return (unsigned int)int64(); }
};

// Faster Uniform Random number generator with period 8.5e37
struct Ranq2 {
  unsigned long long v, w;
  Ranq2(unsigned long long j) : v(4101842887655102017LL), w(1) {
    v ^= j;
    w = int64();
    v = int64();
  }
  inline unsigned long long int64() {
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    return v ^ w;
  }
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  inline unsigned int int32() { return (unsigned int)int64(); }
};

// Implements Knuth’s subtractive generator using only floating operations
struct Ranfib {
  double dtab[55], dd;
  int inext, inextp;

  Ranfib(unsigned long long j) : inext(0), inextp(31) {
    Ranq1 init(j);
    for (int k = 0; k<55; k++) dtab[k] = init.doub();
  }
  // Returns random double-precision floating value between 0. and 1.
  inline double doub() {
    if (++inext == 55) inext = 0;
    if (++inextp == 55) inextp = 0;
    dd = dtab[inext] - dtab[inextp];
    if (dd < 0) dd += 1.0;
    return (dtab[inext] = dd);
  }
  // Returns random 32-bit integer. Recommended only for testing purposes.
  inline unsigned long int32() {
    return (unsigned long)(doub() * 4294967295.0);
  }
};


// Uniform distribution of points on the circumference of a unit circle
// Ref: http://mathworld.wolfram.com/CirclePointPicking.html
template<typename MyRan>
void circle_point_picking(double &x, double &y, MyRan &myran) {
  double a, b, aa, bb, S;
  do {
    a = myran.doub() * 2 - 1;
    b = myran.doub() * 2 - 1;
    aa = a * a;
    bb = b * b;
    S = aa + bb;
  } while (S >= 1);
  x = (aa - bb) / S;
  y = 2 * a * b / S;
}

// Uniform distribution of points on the surface of a unit sphere
// Ref: http://mathworld.wolfram.com/SpherePointPicking.html
template<class MyRan>
void sphere_point_picking(double &x, double &y, double &z, MyRan &myran) {
  double a, b, S;
  do {
    a = myran.doub() * 2 - 1;
    b = myran.doub() * 2 - 1;
    S = a * a + b * b;
  } while (S >= 1);
  double R = std::sqrt(1 - S);
  x = 2 * a * R;
  y = 2 * b * R;
  z = 1 - 2 * S;
}

// Unifor distribution of points on the surface of a unit 4d sphere
// Ref: http://mathworld.wolfram.com/HyperspherePointPicking.html
template<class MyRan>
void hypersphere_point_picking(double *X, MyRan &myran) {
  double a, b, S1, S2;
  do {
    a = myran.doub() * 2 - 1;
    b = myran.doub() * 2 - 1;
    S1 = a * a + b * b;
  } while (S1 >= 1);
  X[0] = a;
  X[1] = b;
  do {
    a = myran.doub() * 2 - 1;
    b = myran.doub() * 2 - 1;
    S2 = a * a + b * b;
  } while (S2 >= 1);
  double Q = std::sqrt((1 - S1) / S2);
  X[2] = a * Q;
  X[3] = b * Q;
}

// For other cases of random point picking, ref.:
// http://mathworld.wolfram.com/topics/RandomPointPicking.html

// Shuffle a array randomly
template<class T, class MyRan>
void shuffle(T *a, int n, MyRan &myran) {
  for (int i = n - 1; i > 0; i--) {
    // generate a random int j that 0 <= j <= i  
    int j = int(myran.doub() * (i + 1));
    if (j > i)
      j = i;
    else if (j < 0)
      j = 0;
    T tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }
}

template<class T, class MyRan>
void shuffle(std::vector<T> &arr, MyRan &myran) {
  const unsigned int n = arr.size();
  for (unsigned int i = n - 1; i > 0; i--) {
    unsigned int j = myran.doub() * (i + 1);
    //auto j = myran.int64() % (i + 1);
    T tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
}

template<class T, class Myran, class UniFunc>
void for_each_shuffle(std::vector<T> &arr, Myran &myran, UniFunc f) {
  const unsigned int n = arr.size();
  for (unsigned int i = n - 1; i > 0; i--) {
    unsigned int j = myran.doub() * (i + 1);
    T tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
    f(arr[i]);
  }
  f(arr[0]);
}
#endif



