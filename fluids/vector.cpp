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




