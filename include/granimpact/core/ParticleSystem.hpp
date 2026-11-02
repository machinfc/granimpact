#pragma once
// Particle storage in SoA (structure of arrays) format.
// Why SoA: it is what lets OpenMP and above all CUDA work later
// without rewriting anything; if this were a vector<struct Particle>, the GPU
// stage would force the whole core to be redone.
#include <vector>

#include "granimpact/core/Types.hpp"

namespace granimpact::core {

class ParticleSystem {
public:
    explicit ParticleSystem(std::size_t capacity = 0) { reserve(capacity); }

    void reserve(std::size_t n);
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

    // Add a particle. Computes mass and inertia from radius and density.
    Index add(const Vec3& position, Real radius, Real density, ParticleKind kind);

    void clear();
    void clearForces();                        // zeroes force, torque, and per-step contact diagnostics
    void applyGravity(Real gz);                // accumulates m*g into each particle's force
    void scaleMassOf(Index i, Real factor);    // rescales mass and inertia of particle i

    // Component accessors (what the kernels consume).
    [[nodiscard]] const std::vector<Vec3>& position() const { return position_; }
    [[nodiscard]] const std::vector<Vec3>& velocity() const { return velocity_; }
    [[nodiscard]] const std::vector<Vec3>& omega() const { return omega_; }
    [[nodiscard]] const std::vector<Vec3>& force() const { return force_; }
    [[nodiscard]] const std::vector<Vec3>& torque() const { return torque_; }
    [[nodiscard]] const std::vector<Real>& radius() const { return radius_; }
    [[nodiscard]] const std::vector<Real>& mass() const { return mass_; }
    [[nodiscard]] const std::vector<Real>& invMass() const { return inv_mass_; }
    [[nodiscard]] const std::vector<Real>& inertia() const { return inertia_; }
    [[nodiscard]] const std::vector<Real>& invInertia() const { return inv_inertia_; }
    [[nodiscard]] const std::vector<ParticleKind>& kind() const { return kind_; }
    [[nodiscard]] const std::vector<Real>& young() const { return young_; }
    [[nodiscard]] const std::vector<Real>& poisson() const { return poisson_; }
    [[nodiscard]] const std::vector<Real>& friction() const { return friction_; }
    [[nodiscard]] const std::vector<Real>& restitution() const { return restitution_; }
    [[nodiscard]] const std::vector<Real>& normalForceSum() const { return normal_force_sum_; }
    [[nodiscard]] const std::vector<Index>& coordination() const { return coordination_; }

    [[nodiscard]] std::vector<Vec3>& position() { return position_; }
    [[nodiscard]] std::vector<Vec3>& velocity() { return velocity_; }
    [[nodiscard]] std::vector<Vec3>& omega() { return omega_; }
    [[nodiscard]] std::vector<Vec3>& force() { return force_; }
    [[nodiscard]] std::vector<Vec3>& torque() { return torque_; }
    [[nodiscard]] std::vector<Real>& radius() { return radius_; }
    [[nodiscard]] std::vector<Real>& invMass() { return inv_mass_; }
    [[nodiscard]] std::vector<Real>& invInertia() { return inv_inertia_; }
    [[nodiscard]] std::vector<ParticleKind>& kind() { return kind_; }
    [[nodiscard]] std::vector<Real>& normalForceSum() { return normal_force_sum_; }
    [[nodiscard]] std::vector<Index>& coordination() { return coordination_; }

    // Global properties (per particle, although uniform today: the interface
    // already supports mixed materials with no changes).
    void setMaterial(Real young, Real poisson, Real friction, Real restitution);
    void setMaterialOf(Index i, Real young, Real poisson, Real friction, Real restitution);

    // Energies and diagnostic magnitudes.
    [[nodiscard]] Real kineticEnergy() const;
    [[nodiscard]] Real rotationalKineticEnergy() const;
    [[nodiscard]] Real gravitationalPotentialEnergy(Real gz) const;
    [[nodiscard]] Real maxRadius() const noexcept { return max_radius_; }
    [[nodiscard]] Real totalMass() const;
    // Solid volume occupied by the grains: the basis of every measured density.
    [[nodiscard]] Real solidsVolume() const;
    [[nodiscard]] Real solidsVolumeOf(ParticleKind k) const;

    // Speed of the fastest particle: used by the stability criterion.
    [[nodiscard]] Real maxSpeed() const;

private:
    std::vector<Vec3> position_, velocity_, omega_, force_, torque_;
    std::vector<Real> radius_, mass_, inv_mass_, inertia_, inv_inertia_;
    std::vector<ParticleKind> kind_;
    std::vector<Real> young_, poisson_, friction_, restitution_;
    // Per-step contact diagnostics (populated by ContactModel::computeForces):
    // - normal_force_sum_[i]: scalar sum of |F_n| acting on particle i [N] (force chains / virial pressure)
    // - coordination_[i]: number of active contacts Z_i on particle i
    std::vector<Real> normal_force_sum_;
    std::vector<Index> coordination_;
    Real max_radius_{0};
    std::size_t size_{0};
};

}  // namespace granimpact::core
