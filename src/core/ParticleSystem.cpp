#include "granimpact/core/ParticleSystem.hpp"

#include <algorithm>
#include <stdexcept>

namespace granimpact::core {

void ParticleSystem::reserve(std::size_t n) {
    position_.reserve(n); velocity_.reserve(n); force_.reserve(n); torque_.reserve(n);
    radius_.reserve(n); mass_.reserve(n); inv_mass_.reserve(n); inertia_.reserve(n);
    kind_.reserve(n); young_.reserve(n); poisson_.reserve(n);
    friction_.reserve(n); restitution_.reserve(n);
}

Index ParticleSystem::add(const Vec3& position, Real radius, Real density, ParticleKind kind) {
    if (radius <= 0) throw std::invalid_argument("ParticleSystem::add: radius must be > 0");
    if (density <= 0) throw std::invalid_argument("ParticleSystem::add: density must be > 0");

    const Real volume = Real(4.0 / 3.0) * kPi * radius * radius * radius;
    const Real mass = density * volume;
    // Solid sphere: I = 2/5 m r^2. The scalar inertia is stored (enough for
    // tangential forces without complex rotation at this stage).
    const Real inertia = Real(0.4) * mass * radius * radius;

    position_.push_back(position);
    velocity_.push_back({});
    force_.push_back({});
    torque_.push_back({});
    radius_.push_back(radius);
    mass_.push_back(mass);
    inv_mass_.push_back(Real(1) / mass);
    inertia_.push_back(inertia);
    kind_.push_back(kind);
    young_.push_back(0); poisson_.push_back(0); friction_.push_back(0); restitution_.push_back(0);
    return static_cast<Index>(size_++);
}

void ParticleSystem::clearForces() {
    std::fill(force_.begin(), force_.end(), Vec3{});
    std::fill(torque_.begin(), torque_.end(), Vec3{});
}

void ParticleSystem::applyGravity(Real gz) {
    if (gz == 0) return;
    for (std::size_t i = 0; i < size_; ++i) {
        force_[i].z += mass_[i] * gz;
    }
}

void ParticleSystem::setMaterial(Real young, Real poisson, Real friction, Real restitution) {
    std::fill(young_.begin(), young_.end(), young);
    std::fill(poisson_.begin(), poisson_.end(), poisson);
    std::fill(friction_.begin(), friction_.end(), friction);
    std::fill(restitution_.begin(), restitution_.end(), restitution);
}

void ParticleSystem::setMaterialOf(Index i, Real young, Real poisson, Real friction, Real restitution) {
    const auto k = static_cast<std::size_t>(i);
    young_[k] = young; poisson_[k] = poisson; friction_[k] = friction; restitution_[k] = restitution;
}

Real ParticleSystem::kineticEnergy() const {
    Real ke = 0;
    for (std::size_t i = 0; i < size_; ++i) ke += Real(0.5) * mass_[i] * velocity_[i].norm2();
    return ke;
}

Real ParticleSystem::gravitationalPotentialEnergy(Real gz) const {
    Real pe = 0;
    for (std::size_t i = 0; i < size_; ++i) pe += mass_[i] * (-gz) * position_[i].z;
    return pe;
}

Real ParticleSystem::maxRadius() const {
    Real r = 0;
    for (std::size_t i = 0; i < size_; ++i) r = std::max(r, radius_[i]);
    return r;
}

Real ParticleSystem::totalMass() const {
    Real m = 0;
    for (std::size_t i = 0; i < size_; ++i) m += mass_[i];
    return m;
}

Real ParticleSystem::solidsVolume() const {
    Real volume = 0;
    for (std::size_t i = 0; i < size_; ++i) {
        volume += Real(4.0 / 3.0) * kPi * radius_[i] * radius_[i] * radius_[i];
    }
    return volume;
}

Real ParticleSystem::maxSpeed() const {
    Real v = 0;
    for (std::size_t i = 0; i < size_; ++i) v = std::max(v, velocity_[i].norm());
    return v;
}

}  // namespace granimpact::core
