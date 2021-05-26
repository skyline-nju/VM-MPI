#ifndef VECT_H
#define VECT_H
#include <iostream>
#include <cmath>
#include "rand.h"

template <typename T>
struct Vec_2 {
  T x;
  T y;

  Vec_2() = default;
  Vec_2(T x0, T y0) : x(x0), y(y0) {}
  Vec_2(const Vec_2<T> &a) : x(a.x), y(a.y) {}

  T operator[](int i) const { return *(&x + i); }
  T& operator [](int i) { return *(&x + i); }

  Vec_2<T>& operator =(const Vec_2<T>& rhs);
  Vec_2<T> operator -() const;
  void operator +=(const Vec_2<T>& rhs);
  void operator -=(const Vec_2<T>& rhs);
  void operator *=(double rhs);
  void operator /=(double rhs);
  Vec_2<T> operator +(const Vec_2<T>& rhs) const;
  Vec_2<T> operator -(const Vec_2<T>& rhs) const;
  Vec_2<T> operator +(double rhs) const;
  Vec_2<T> operator -(double rhs) const;
  Vec_2<T> operator *(double rhs) const;
  Vec_2<T> operator /(double rhs) const;

  template <typename T2>
  Vec_2<T> operator * (const Vec_2<T2>& rhs) const;

  double dot(const Vec_2<T>& a) const;
  double cross(const Vec_2<T>& a) const;
  void normalize();
  double square() const;
  double module() const { return std::sqrt(square()); }
  Vec_2 inverse() const;
  void rotate(double theta);

  friend std::ostream & operator << (std::ostream &out, const Vec_2<T> & obj) {
    out << obj.x << "\t" << obj.y;
    return out;
  }
};

template <typename T>
Vec_2<T>& Vec_2<T>::operator=(const Vec_2<T>& rhs) {
  this->x = rhs.x;
  this->y = rhs.y;
  return *this;
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-() const {
  return Vec_2<T>(-x, -y);
  
}

template <typename T>
void Vec_2<T>::operator+=(const Vec_2<T>& rhs) {
  x += rhs.x;
  y += rhs.y;
}

template <typename T>
void Vec_2<T>::operator-=(const Vec_2<T>& rhs) {
  x -= rhs.x;
  y -= rhs.y;
}

template <typename T>
void Vec_2<T>::operator*=(double rhs) {
  x *= rhs;
  y *= rhs;
}

template <typename T>
void Vec_2<T>::operator/=(double rhs) {
  x /= rhs;
  y /= rhs;
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(double rhs) const {
  return Vec_2<T>(x + rhs, y + rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x + rhs.x, y + rhs.y);
}

template <typename T>
Vec_2<T> operator+(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs + rhs.x, lhs + rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(double rhs) const {
  return Vec_2<T>(x - rhs, y - rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x - rhs.x, y - rhs.y);
}

template <typename T>
Vec_2<T> operator-(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs - rhs.x, lhs - rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator*(double rhs) const {
  return Vec_2<T>(x * rhs, y * rhs);
}

template <typename T>
Vec_2<T> operator*(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs * rhs.x, lhs * rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator/(double rhs) const {
  return Vec_2<T>(x / rhs, y / rhs);
}

template <typename T>
template <typename T2>
Vec_2<T> Vec_2<T>::operator * (const Vec_2<T2>& rhs) const {
  return Vec_2<T>(x * rhs.x, y * rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::inverse() const {
  return Vec_2<T>(1 / x, 1 / y);
}

template <typename T>
double Vec_2<T>::square() const {
  return x * x + y * y;
}

template <typename T>
double Vec_2<T>::dot(const Vec_2<T>& a) const {
  return x * a.x + y * a.y;
}

template <typename T>
double Vec_2<T>::cross(const Vec_2<T>& a) const {
  return x * a.y - y * a.x;
}

template <typename T>
void Vec_2<T>::normalize() {
  const auto one_over_r = 1 / std::sqrt(square());
  x *= one_over_r;
  y *= one_over_r;
}

template <typename T>
void Vec_2<T>::rotate(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  const auto x_new = x * c - y * s;
  const auto y_new = x * s + y * c;
  x = x_new;
  y = y_new;
}

template <typename T>
struct Vec_3 {
  T x;
  T y;
  T z;

  Vec_3() = default;
  Vec_3(T x0, T y0, T z0) : x(x0), y(y0), z(z0) {}
  Vec_3(const Vec_3<T> &a) : x(a.x), y(a.y), z(a.z) {}

  T operator[](int i) const { return *(&x + i); }
  T& operator [](int i) { return *(&x + i); }

  Vec_3<T>& operator =(const Vec_3<T>& rhs);
  Vec_3<T> operator -()const;
  void operator +=(const Vec_3<T>& rhs);
  void operator -=(const Vec_3<T>& rhs);
  void operator *=(double rhs);
  void operator /=(double rhs);
  Vec_3<T> operator +(const Vec_3<T>& rhs) const;
  Vec_3<T> operator -(const Vec_3<T>& rhs) const;
  Vec_3<T> operator +(double rhs) const;
  Vec_3<T> operator -(double rhs) const;
  Vec_3<T> operator *(double rhs) const;
  Vec_3<T> operator /(double rhs) const;

  template <typename T2>
  Vec_3<T> operator *(const Vec_3<T2> & rhs) const;

  double dot(const Vec_3<T>& a) const;
  Vec_3<T> cross(const Vec_3<T>& a) const;
  void normalize();
  double square() const;
  double module() const;
  Vec_3 inverse() const;
  void rotate(double theta);
  void rotate(double theta, const Vec_3<double> &a);
  void rotate(double c, double s, const Vec_3<double> &a);

  template <typename TRan>
  void get_perp_vec(Vec_3<T> &v_perp, TRan &myran) const;

  template <typename TRan>
  void rotate_rand(double theta, TRan &myran);


  friend std::ostream & operator << (std::ostream &out, const Vec_3<T> & obj) {
    out << obj.x << "\t" << obj.y << "\t" << obj.z;
    return out;
  }
};

template <typename T>
Vec_3<T>& Vec_3<T>::operator=(const Vec_3<T>& rhs) {
  this->x = rhs.x;
  this->y = rhs.y;
  this->z = rhs.z;
  return *this;
}

template <typename T>
Vec_3<T> Vec_3<T>::operator-() const {
  return Vec_3<T>(-x, -y, -z);
}

template <typename T>
void Vec_3<T>::operator+=(const Vec_3<T>& rhs) {
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
}

template <typename T>
void Vec_3<T>::operator-=(const Vec_3<T>& rhs) {
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
}

template <typename T>
void Vec_3<T>::operator*=(double rhs) {
  x *= rhs;
  y *= rhs;
  z *= rhs;
}

template <typename T>
void Vec_3<T>::operator/=(double rhs) {
  x /= rhs;
  y /= rhs;
  z /= rhs;
}

template <typename T>
Vec_3<T> Vec_3<T>::operator+(double rhs) const {
  return Vec_3<T>(x + rhs, y + rhs, z + rhs);
}

template <typename T>
Vec_3<T> Vec_3<T>::operator+(const Vec_3<T>& rhs) const {
  return Vec_3<T>(x + rhs.x, y + rhs.y, z + rhs.z);
}

template <typename T>
Vec_3<T> operator+(double lhs, const Vec_3<T>& rhs) {
  return Vec_3<T>(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z);
}

template <typename T>
Vec_3<T> Vec_3<T>::operator-(double rhs) const {
  return Vec_3<T>(x - rhs, y - rhs, z - rhs);
}

template <typename T>
Vec_3<T> Vec_3<T>::operator-(const Vec_3<T>& rhs) const {
  return Vec_3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
}

template <typename T>
Vec_3<T> operator-(double lhs, const Vec_2<T>& rhs) {
  return Vec_3<T>(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
}

template <typename T>
Vec_3<T> Vec_3<T>::operator*(double rhs) const {
  return Vec_3<T>(x * rhs, y * rhs, z * rhs);
}

template <typename T>
Vec_3<T> operator*(double lhs, const Vec_3<T>& rhs) {
  return Vec_3<T>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
}

template <typename T>
Vec_3<T> Vec_3<T>::operator/(double rhs) const {
  return Vec_3<T>(x / rhs, y / rhs, z / rhs);
}

template <typename T>
template <typename T2>
Vec_3<T> Vec_3<T>::operator*(const Vec_3<T2>& rhs) const {
  return Vec_3<T>(x * rhs.x, y * rhs.y, z * rhs.z);
}

template <typename T>
Vec_3<T> Vec_3<T>::inverse() const {
  return Vec_3<T>(1 / x, 1 / y, 1 / z);
}

template <typename T>
double Vec_3<T>::square() const {
  return x * x + y * y + z * z;
}

template <typename T>
double Vec_3<T>::module() const {
  return std::sqrt(square());
}

template <typename T>
double Vec_3<T>::dot(const Vec_3<T>& a) const {
  return x * a.x + y * a.y + z * a.z;
}

template <typename T>
Vec_3<T> Vec_3<T>::cross(const Vec_3<T>& a) const {
  return Vec_3<T>(y * a.z - z * a.y,
                  z * a.x - x * a.z,
                  x * a.y - y * a.x);
}

template <typename T>
void Vec_3<T>::normalize() {
  const auto one_over_r = 1 / std::sqrt(square());
  x *= one_over_r;
  y *= one_over_r;
  z *= one_over_r;
}

/**
 * @brief Rotate around the z axis (0, 0, 1) by theta
 * 
 * @tparam T      
 * @param theta   The angle need to rotate
 */
template <typename T>
void Vec_3<T>::rotate(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  const auto x_new = x * c - y * s;
  const auto y_new = x * s + y * c;
  x = x_new;
  y = y_new;
}

/**
 * @brief Rotate around a vector by a given angle
 * 
 * @tparam T 
 * @param theta  Thg angle to rotate
 * @param a      Rotation axis
 */
template <typename T>
void Vec_3<T>::rotate(double theta, const Vec_3<double> &a) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  rotate(c, s, a);
}

/**
 * @brief    Rotate aound a vector by a given angle
 * 
 * @tparam T 
 * @param c   cos theta
 * @param s   sin theta
 * @param a   Rotaton axis, should be normalized.
 */
template <typename T>
void Vec_3<T>::rotate(double c, double s, const Vec_3<double> &a) {
  const auto bxx = a.x * a.x * (1 - c);
  const auto bxy = a.x * a.y * (1 - c);
  const auto bxz = a.x * a.z * (1 - c);
  const auto byy = a.y * a.y * (1 - c);
  const auto byz = a.y * a.z * (1 - c); 
  const auto bzz = a.z * a.z * (1 - c);
  const auto x_new = (bxx + c) * x + (bxy - a.z * s) * y + (bxz + a.y * s) * z;
  const auto y_new = (bxy + a.z * s) * x + (byy + c) * y + (byz - a.x * s) * z;
  const auto z_new = (bxz - a.y * s) * x + (byz + a.x * s) * y + (bzz + c) * z;
  x = x_new;
  y = y_new;
  z = z_new;
}

template <typename T>
template <typename TRan>
void Vec_3<T>::get_perp_vec(Vec_3<T>& v_perp, TRan& myran) const {
  double s = std::sqrt(x * x + y * y);
  Vec_3<T> rot_axis(-y / s, x / s, 0.);
  double c = z;

  circle_point_picking(v_perp.x, v_perp.y, myran);
  v_perp.z = 0;
  v_perp.rotate(c, s, rot_axis);
}

/**
 * @brief Rotate around a random axis by angle theta
 * 
 * The vector to be rotated should be normalized.
 * 
 * @tparam T     
 * @tparam TRan  Template for random number generator.
 * @param theta  Angle to rotate.
 * @param myran  Random number generator
 */
template <typename T>
template <typename TRan>
void Vec_3<T>::rotate_rand(double theta, TRan& myran) {
  Vec_3<double> v_perp{};
  get_perp_vec(v_perp, myran);
  rotate(theta, v_perp);
}

#endif


