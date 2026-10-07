#pragma once
// Base types of the simulator. All the code uses `Real` so it can switch to
// float (or to mixed precision on GPU) without touching the physics.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace granimpact::core {

using Real = double;
using Index = std::int32_t;
inline constexpr Real kPi = 3.14159265358979323846;

struct Vec3 {
    Real x{}, y{}, z{};

    constexpr Vec3() = default;
    constexpr Vec3(Real x_, Real y_, Real z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator*(Real s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(Real s) const { return {x / s, y / s, z / s}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(Real s) { x *= s; y *= s; z *= s; return *this; }

    [[nodiscard]] constexpr Real dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    [[nodiscard]] constexpr Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    [[nodiscard]] Real norm() const { return std::sqrt(dot(*this)); }
    [[nodiscard]] Real norm2() const { return dot(*this); }
    [[nodiscard]] Vec3 normalized() const {
        const Real n = norm();
        return n > Real(0) ? (*this) / n : Vec3{};
    }
};

struct Box {
    Vec3 min{};
    Vec3 max{};

    [[nodiscard]] Vec3 size() const { return max - min; }
    [[nodiscard]] Real volume() const { return size().x * size().y * size().z; }
    [[nodiscard]] Vec3 center() const { return (min + max) * Real(0.5); }
};

enum class ParticleKind : std::uint8_t { Bed = 0, Projectile = 1, Ejecta = 2 };

}  // namespace granimpact::core
