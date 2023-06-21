#include "vector.h"
#include <cmath>
#include <stdexcept>
#include <ostream>

// NOTE: No need to call std::move anywhere. 
// See this: https://stackoverflow.com/questions/14856344/when-should-stdmove-be-used-on-a-function-return-value

Vector3 operator+(const Vector3& a, const Vector3& b) {
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector3 operator-(const Vector3& a, const Vector3& b) {
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3 operator*(const double a, const Vector3& b) {
    return Vector3(a*b.x, a*b.y, a*b.z);
}

Vector3 operator*(const Vector3& a, const double b) {
    return Vector3(a.x*b, a.y*b, a.z*b);
}

Vector3 operator/(const Vector3& a, const double b) {
    return Vector3(a.x / b, a.y / b, a.z / b);
}

double Vector3::norm2() const {
    return x * x + y * y + z * z;
}

double Vector3::norm() const {
    return sqrt(norm2());
}

double Vector3::dot(const Vector3& b) const {
    return x * b.x + y * b.y + z * b.z;
}

Vector3 Vector3::cross(const Vector3& b) const {
    return Vector3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
}

Vector3 Vector3::operator-() const {
    return Vector3(-x, -y, -z);
}

inline double gamma_correct(double x) {
    return std::pow(x, 1.0 / 2.2);
}

Vector3 gamma_correct(const Vector3& v) {
    return Vector3(gamma_correct(v.x), gamma_correct(v.y), gamma_correct(v.z));
}

inline double reverse_gamma_correct(double x) {
    return std::pow(x, 2.2);
}

Vector3 reverse_gamma_correct(const Vector3& v) {
    return Vector3(reverse_gamma_correct(v.x), reverse_gamma_correct(v.y), reverse_gamma_correct(v.z));
}

Vector3 Vector3::normalized() const {
    if (norm2() <= 0.0001)
        return Vector3(0, 0, 0);
    else 
        return *this / norm();
}

double Vector3::dist2(const Vector3& b) const {
    return (*this - b).norm2();
}

double Vector3::dist(const Vector3& b) const {
    return (*this - b).norm();
}

Vector3 operator*(const Vector3& a, const Vector3& b) {
    return Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
}

void Vector3::normalize() {
    double n = norm();
    if (n > 0.0001) {
        x /= n;
        y /= n;
        z /= n;
    }
}

void Vector3::operator+=(const Vector3& b) {
    x += b.x;
    y += b.y;
    z += b.z;
}

double& Vector3::operator[](size_t i) {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    throw std::out_of_range("Vector3 index out of range");
}

std::ostream& operator<< ( std::ostream& os, const Vector3& c ) {
    os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
    return os;
}