#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <array>
#include <ostream>

class Vector2 {
public:
    double x, y;
    Vector2(double x = 0, double y = 0): x(x), y(y) {};
    Vector2(const std::array<double, 2>& a): x(a[0]), y(a[1]) {};

	double norm2() const;
	double norm() const;

    double dot(const Vector2& b) const;
    double dist2(const Vector2& b) const;
    double dist(const Vector2& b) const;
    double cross(const Vector2& b) const;

    Vector2 normalized() const;
    Vector2 operator-() const;

    void normalize();
    void operator+=(const Vector2& b);
    double& operator[](size_t i);
};

Vector2 operator+(const Vector2& a, const Vector2& b);
Vector2 operator-(const Vector2& a, const Vector2& b);
Vector2 operator*(const double a, const Vector2& b);
Vector2 operator*(const Vector2& a, const double b);
Vector2 operator/(const Vector2& a, const double b);
Vector2 operator*(const Vector2& a, const Vector2& b);

Vector2 gamma_correct(const Vector2& v);
Vector2 reverse_gamma_correct(const Vector2& v);

#endif