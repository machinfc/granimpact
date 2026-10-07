#include "granimpact/sim/Simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include "granimpact/generation/Generators.hpp"
#include "granimpact/physics/Integrator.hpp"

namespace granimpact::sim {

using core::ParticleKind;
using core::Real;
using core::SimConfig;
using core::Vec3;

namespace {

// Initial projectile kinetic energy: it is the E that enters the models and the
// one reported in the CSV. It is measured, not estimated.
Real projectileKineticEnergy(const core::ParticleSystem& particles) {
    Real energy = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Projectile) continue;
        energy += Real(0.5) * particles.mass()[i] * particles.velocity()[i].norm2();
    }
    return energy;
}

// z bound of the projectile centre of mass (0 if there is no projectile).
Real projectileComZ(const core::ParticleSystem& particles) {
    Real mass = 0, mz = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Projectile) continue;
        mass += particles.mass()[i];
        mz += particles.position()[i].z * particles.mass()[i];
    }
    return mass > 0 ? mz / mass : 0;
}

// Velocity of the projectile centre of mass: this is the REAL impact velocity,
// which may differ from the config one because of the free fall before contact.
Real projectileComSpeed(const core::ParticleSystem& particles) {
    Real mass = 0;
    Vec3 momentum{};
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Projectile) continue;
        mass += particles.mass()[i];
        momentum += particles.velocity()[i] * particles.mass()[i];
    }
    return mass > 0 ? (momentum / mass).norm() : 0;
}

}  // namespace

Simulation::Simulation(SimConfig config, std::string backend)
    : config_(std::move(config)),
      solver_(backends::makeContactSolver(backend, config_)),   // same order as the class
      analyzer_(config_) {}

std::pair<Real, bool> Simulation::rayleighStabilityCheck() const {
    // Rayleigh time for Hertz spheres (Zhang & Makse):
    //   t_R = pi * r * sqrt(rho / G) / (0.1631 * nu + 0.8766)
    // dt <= t_R / 20 is required so the contact is not integrated "in jumps".
    const auto& pc = config_.particles;
    const Real g = pc.young / (2.0 * (1.0 + pc.poisson));
    const Real t_r = core::kPi * pc.r_max * std::sqrt(pc.rho_grain / g)
                     / (0.1631 * pc.poisson + 0.8766);
    return {t_r / 20.0, config_.simulation.dt <= t_r / 20.0};
}

io::StepRecord Simulation::step(io::StepRecord record, bool measure_crater, Real gravity_scale) {
    physics::VelocityVerlet integrator(particles_, config_.simulation.dt, config_.gravity.gz);
    const core::Box box = config_.domain.box();
    const Real gz = config_.gravity.gz * gravity_scale;

    // 1) kick / drift / fronteras
    integrator.firstHalfKick(particles_);
    integrator.drift(particles_);
    physics::applyBoundaries(particles_, box, config_);

    // 2) forces at the new position
    particles_.clearForces();
    particles_.applyGravity(gz);
    const Real cell_size = physics::CellList::suggestCellSize(particles_.maxRadius());
    cells_.build(particles_, box, cell_size);
    const auto stats = solver_->computeForces(particles_, cells_, config_.simulation.dt);

    // 3) second half of the kick
    integrator.secondHalfKick(particles_);

    // Impact velocity: measured when the projectile centre of mass reaches the
    // bed surface (geometric criterion). Using contacts does not work, because
    // during settling and the impact phase itself there are grain-grain contacts
    // that have nothing to do with the projectile.
    if (!first_contact_seen_ && analyzer_.hasReference() && z_surface_ > 0) {
        const Real com = projectileComZ(particles_);
        if (com > 0 && com - config_.projectile.radius
                           <= z_surface_ + 2 * config_.particles.r_max) {
            first_contact_seen_ = true;
            impact_speed_ = projectileComSpeed(particles_);
        }
    }

    record.time += config_.simulation.dt;
    record.particles = particles_.size();
    record.contacts = stats.touching;
    record.kinetic_energy = particles_.kineticEnergy();
    record.potential_energy = particles_.gravitationalPotentialEnergy(gz);
    record.max_speed = particles_.maxSpeed();
    if (measure_crater) {
        record.crater = analyzer_.analyze(particles_, impact_energy_);
        record.has_crater = true;
    }
    return record;
}

SimulationResult Simulation::run(std::ostream& log) {
    const auto wall_start = std::chrono::steady_clock::now();
    const auto [dt_max, dt_ok] = rayleighStabilityCheck();
    log << "GranImpact - " << config_.simulation.name << "  (backend " << solver_->name() << ")\n";
    log << "  dt = " << config_.simulation.dt << " s  (dt_max Rayleigh " << dt_max << " s)"
        << (dt_ok ? " [stable]" : " [WARNING: dt above the limit, it may blow up]") << "\n";

    std::mt19937_64 rng(config_.simulation.seed);

    // ---- PHASE 1: granular bed ------------------------------------------------
    const auto bed = generation::generateBed(particles_, config_, rng);
    log << "  Bed: " << bed.count << " grains seeded in " << bed.region_height
        << " m  (lattice " << bed.lattice_spacing * 1000 << " x " << bed.layer_spacing * 1000
        << " mm, seeded phi " << bed.seeded_phi << ", target " << bed.target_phi << ")\n";

    const auto out_dir = std::filesystem::path(config_.output.dir);
    std::filesystem::create_directories(out_dir);
    if (config_.output.csv) {
        csv_ = std::make_unique<io::CsvWriter>(
            out_dir / (config_.simulation.name + "_crater_evolution.csv"), config_);
    }
    if (config_.output.vtk) {
        vtk_ = std::make_unique<io::VtkWriter>(out_dir, config_.simulation.name);
    }
    if (vtk_) vtk_->write(particles_, 0, 0.0);

    // ---- PHASE 2: settling ----------------------------------------------------
    std::uint64_t settle_steps = static_cast<std::uint64_t>(
        std::llround(config_.simulation.t_settle / config_.simulation.dt));
    if (config_.settle.max_steps > 0) settle_steps = config_.settle.max_steps;
    log << " Settling: up to " << settle_steps << " steps"
        << (config_.settle.gravity_boost != 1.0
                ? "  (gravedad x" + std::to_string(config_.settle.gravity_boost) + ")"
                : "")
        << " (cut by energy convergence)\n";

    io::StepRecord record;
    const std::uint64_t check_every = std::max<std::uint64_t>(settle_steps / 200, 1);
    std::uint64_t used_steps = 0;
    bool converged = false;
    // Convergence criterion on the bed HEIGHT, not only on kinetic energy:
    // a bed can have Ek ~ 0 and still creep (filling gaps of the initial
    // packing). Cutting there makes that subsidence show up later as a false
    // "crater". The criterion requires the mean height to stop changing.
    const Real damping = config_.settle.velocity_damping;
    const auto meanBedHeight = [&]() {
        Real sum = 0;
        std::size_t n = 0;
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            if (particles_.kind()[i] != core::ParticleKind::Bed) continue;
            sum += particles_.position()[i].z;
            ++n;
        }
        return n > 0 ? sum / static_cast<Real>(n) : 0;
    };

    Real height_prev = particles_.empty() ? 0 : meanBedHeight();
    const Real height_tolerance = Real(1.0e-3) * config_.particles.r_mean();
    for (std::uint64_t s = 0; s < settle_steps; ++s) {
        record = step(record, false, config_.settle.gravity_boost);
        if (damping > 0) {
            const Real factor = std::exp(-damping * config_.simulation.dt);
            for (core::Vec3& v : particles_.velocity()) v = v * factor;
        }
        used_steps = s + 1;
        if (s + 1 >= config_.settle.min_steps && s % check_every == 0 && !particles_.empty()) {
            const Real ke_per_particle = particles_.kineticEnergy()
                                         / static_cast<Real>(particles_.size());
            const Real height_now = meanBedHeight();
            const bool height_stable = std::abs(height_now - height_prev) < height_tolerance;
            height_prev = height_now;
            if (ke_per_particle < config_.settle.ke_tolerance && height_stable) {
                converged = true;
                break;
            }
        }
    }
    settle_steps_used_ = used_steps;
    settle_converged_ = converged;
    log << " Settling finished in " << used_steps << " steps (t = " << record.time
        << " s, Ek/particula = "
        << (particles_.empty() ? 0.0
                                : particles_.kineticEnergy() / static_cast<core::Real>(particles_.size()))
           << " J)"
        << (converged ? "  [convergido]" : "  [NO convergio: sube t_settle]") << "\n";

    // Reference surface = the already settled bed: it is the crater's zero.
    analyzer_.recordReferenceSurface(particles_);
    z_surface_ = analyzer_.referenceHeight();
    const Real settled_phi = particles_.solidsVolume()
                             / (config_.domain.x * config_.domain.y
                                * std::max(analyzer_.referenceHeight(), Real(1e-9)));
    log << " Settled bed: z_ref = " << analyzer_.referenceHeight() << " m, measured phi "
        << settled_phi << "\n";

    // ---- PHASE 3: projectile --------------------------------------------------
    // Safe placement: the projectile cannot be born inside the bed. If the height
    // in the config leaves no clearance, it is raised and a warning is issued (they never overlap).
    const Real z_surface = z_surface_;
    const Real z_min_allowed = z_surface + config_.projectile.radius + 2 * config_.particles.r_max;
    if (config_.projectile.start_z < z_min_allowed) {
        log << "  WARNING: start_z = " << config_.projectile.start_z
            << " m falls inside the bed (surface " << z_surface << " m). Raising to "
            << z_min_allowed << " m.\n";
        config_.projectile.start_z = z_min_allowed;
    }
    const auto proj = generation::generateProjectile(particles_, config_, rng);
    impact_energy_ = projectileKineticEnergy(particles_);
    // Free fall before contact: estimated to warn if t_end falls short.
    const Real gap = config_.projectile.start_z - config_.projectile.radius - z_surface;
    const Real v0 = config_.projectile.speed;
    const Real fall_time = gap > 0 && std::abs(config_.gravity.gz) > 0
                               ? (std::sqrt(v0 * v0 + 2 * std::abs(config_.gravity.gz) * gap) - v0)
                                     / std::abs(config_.gravity.gz)
                               : 0;
    log << " Projectile: " << proj << " grains at " << config_.projectile.speed << " m/s and "
        << config_.projectile.angle_deg << " degrees | launch E = " << impact_energy_
        << " J\n";
    log << " Gap to the bed " << gap << " m -> impact at ~" << fall_time << " s (t_end "
        << config_.simulation.t_end << " s)"
        << (fall_time + 0.02 > config_.simulation.t_end ? "  [WARNING: raise t_end]" : "") << "\n";

    if (vtk_) vtk_->write(particles_, 1, record.time);

    // ---- PHASE 4: impact -------------------------------------------------------
    const auto impact_steps = static_cast<std::uint64_t>(
        std::llround(config_.simulation.t_end / config_.simulation.dt));
    log << " Impact: " << impact_steps << " steps...\n";

    std::uint64_t vtk_step = 2;
    const auto progress_every = std::max<std::uint64_t>(impact_steps / 10, 1);
    for (std::uint64_t s = 0; s < impact_steps; ++s) {
        record = step(record, false);

        if (csv_ && (s % config_.simulation.output_every == 0)) {
            record.crater = analyzer_.analyze(particles_, impact_energy_);
            record.has_crater = true;
            csv_->writeRow(record);
        }
        if (vtk_ && vtk_step % config_.output.vtk_every == 0) {
            vtk_->write(particles_, vtk_step, record.time);
        }
        ++vtk_step;

        if (s % progress_every == 0) {
            const auto c = analyzer_.analyze(particles_, impact_energy_);
            log << "    t=" << std::fixed << std::setprecision(4) << record.time << " s | D="
                << std::setprecision(1) << c.D * 1000.0 << " mm | d_exc=" << c.d_exc * 1000.0
                << " mm | h_rim=" << c.h_rim * 1000.0 << " mm | Z/D=" << std::setprecision(3)
                << c.aspect_zd << " | eps=" << c.epsilon << " | " << c.classification << "\n";
        }
    }

    // ---- Final result ----------------------------------------------------------
    SimulationResult result;
    result.crater = analyzer_.analyze(particles_, impact_energy_);
    result.last_step = record;
    result.simulated_time = record.time;
    result.bed_particles = bed.count;
    result.projectile_particles = proj;
    result.impact_energy = impact_energy_;
    result.impact_speed = impact_speed_;
    result.settle_steps = settle_steps_used_;
    result.settle_converged = settle_converged_;
    result.seeded_phi = bed.seeded_phi;
    result.settled_phi = settled_phi;

    if (csv_) {
        record.crater = result.crater;
        record.has_crater = true;
        csv_->writeRow(record);
        result.csv_path = csv_->path();
    }
    {
        io::CaseSummary summary;
        summary.backend = solver_->name();
        summary.bed_particles = bed.count;
        summary.projectile_particles = proj;
        summary.settle_steps = settle_steps_used_;
        summary.settle_converged = settle_converged_;
        summary.impact_speed = impact_speed_;
        summary.impact_energy = impact_energy_;
        summary.seeded_phi = bed.seeded_phi;
        summary.settled_phi = settled_phi;
        summary.simulated_time = record.time;
        summary.wall_seconds =
            std::chrono::duration<core::Real>(std::chrono::steady_clock::now() - wall_start).count();
        const auto path = std::filesystem::path(config_.output.dir)
                           / (config_.simulation.name + "_summary.txt");
        io::writeCaseSummary(path, config_, result.crater, summary);
        result.summary_path = path;
        if (csv_) csv_->writeFooter(config_, result.crater, summary);
    }
    if (vtk_) {
        vtk_->write(particles_, vtk_step, record.time);
        vtk_->writeCollection(config_);
        result.vtk_collection = std::filesystem::path(config_.output.dir)
                                / (config_.simulation.name + ".pvd");
    }

    log << "\n  Final result\n";
    log << "    D         = " << result.crater.D * 1000.0 << " mm\n";
    log << "    d_exc     = " << result.crater.d_exc * 1000.0 << " mm\n";
    log << "    h_rim     = " << result.crater.h_rim * 1000.0 << " mm\n";
    log << "    epsilon   = " << result.crater.epsilon << "\n";
    log << "    Z/D       = " << result.crater.aspect_zd << "\n";
    log << "    V_in      = " << result.crater.v_in * 1e6 << " cm3\n";
    log << "    Class     = " << result.crater.classification << "\n";
    log << " Bed = phi " << result.crater.phi_bed << " measured, height "
        << result.crater.bed_height * 1000.0 << " mm\n";
    log << " v impact  = " << result.impact_speed << " m/s (measured at first contact)\n";
    log << "    E impact  = " << result.crater.impact_energy << " J | ejecta "
        << result.crater.ejecta_mass * 1000.0 << " g at "
        << result.crater.ejecta_speed_mean << " m/s mean\n";
    log << "    Models    : D_log = " << result.crater.d_log_model * 1000.0
        << " mm | D_uehara = " << result.crater.d_uehara_model * 1000.0
        << " mm | d_heckel = " << result.crater.d_heckel_model * 1000.0 << " mm\n";
    return result;
}

// z bound of the projectile centre of mass (0 if there is no projectile).
Real projectileComZ(const core::ParticleSystem& particles) {
    Real mass = 0, mz = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Projectile) continue;
        mass += particles.mass()[i];
        mz += particles.position()[i].z * particles.mass()[i];
    }
    return mass > 0 ? mz / mass : 0;
}

// Velocity of the projectile centre of mass: this is the REAL impact velocity,
// which may differ from the config one because of the free fall before contact.
Real projectileComSpeed(const core::ParticleSystem& particles) {
    Real mass = 0;
    Vec3 momentum{};
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Projectile) continue;
        mass += particles.mass()[i];
        momentum += particles.velocity()[i] * particles.mass()[i];
    }
    return mass > 0 ? (momentum / mass).norm() : 0;
}

}  // namespace granimpact::sim
