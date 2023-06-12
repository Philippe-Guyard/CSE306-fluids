#include "vector.h"
#include <cmath>
#include <stdexcept>
#include <ostream>

double Vector2::norm2() const {
    return x*x + y*y;
}

Vector2 operator+(const Vector2& a, const Vector2& b) {
    return Vector2(a.x + b.x, a.y + b.y);
}

Vector2 operator-(const Vector2& a, const Vector2& b) {
    return Vector2(a.x - b.x, a.y - b.y);
}

Vector2 operator*(double a, const Vector2& b) {
    return Vector2(a * b.x, a * b.y);
}

Vector2 operator*(const Vector2& a, double b) {
    return Vector2(a.x * b, a.y * b);
}

Vector2 operator/(const Vector2& a, double b) {
    return Vector2(a.x / b, a.y / b);
}

double Vector2::dot(const Vector2& other) const {
    return x * other.x + y * other.y;
}   

double Vector2::cross(const Vector2& other) const {
    return x * other.y - y * other.x;
}

double Vector2::norm() const {
    return std::sqrt(norm2());
}

Vector2 Vector2::normalized() const {
    return *this / norm();
}

std::ostream& operator<<(std::ostream& os, const Vector2& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

void Vector2::operator+=(const Vector2& b) {
    x += b.x;
    y += b.y;
}

double& Vector2::operator[](size_t i) {
    if(i == 0) return x;
    if(i == 1) return y;
    throw std::runtime_error("Vector2 has only 2 elements");
}

Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}

double Vector2::dist2(const Vector2& b) const {
    return (*this - b).norm2();
}

double Vector2::dist(const Vector2& b) const {
    return (*this - b).norm();
}

Vector2 operator*(const Vector2& a, const Vector2& b) {
    return Vector2(a.x * b.x, a.y * b.y);
}

Vector2 gamma_correct(const Vector2& v) {
    return Vector2(std::sqrt(v.x), std::sqrt(v.y));
}

Vector2 reverse_gamma_correct(const Vector2& v) {
    return Vector2(v.x * v.x, v.y * v.y);
}

void Vector2::normalize() {
    x /= norm();
    y /= norm();
}



