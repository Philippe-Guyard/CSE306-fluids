#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <array>
#include <ostream>

class Vector3 {
public:
    double x, y, z;
    Vector3(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z) {};
    // Not necessary since C++ has a default copy constructor
    // Vector3(const Vector3& v): x(v.x), y(v.y), z(v.z) {};
    Vector3(const std::array<double, 3>& a): x(a[0]), y(a[1]), z(a[2]) {};

	double norm2() const;
	double norm() const;

    double dot(const Vector3& b) const;
    double dist2(const Vector3& b) const;
    double dist(const Vector3& b) const;
    Vector3 cross(const Vector3& b) const;

    Vector3 normalized() const;
    Vector3 operator-() const;

    void normalize();
    void operator+=(const Vector3& b);
    double& operator[](size_t i);
};

Vector3 operator+(const Vector3& a, const Vector3& b);
Vector3 operator-(const Vector3& a, const Vector3& b);
Vector3 operator*(const double a, const Vector3& b);
Vector3 operator*(const Vector3& a, const double b);
Vector3 operator/(const Vector3& a, const double b);
Vector3 operator*(const Vector3& a, const Vector3& b);

Vector3 gamma_correct(const Vector3& v);
Vector3 reverse_gamma_correct(const Vector3& v);

std::ostream& operator<< ( std::ostream& os, const Vector3& c );

#endif