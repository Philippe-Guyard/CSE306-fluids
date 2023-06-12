#pragma once
#include <vector>

#include "vector.h"
#include "ot.hpp"
#include "poly.hpp"

#define G 9.81

const double DEFAULT_MASS_PARTICLE = 200.;
const double DEFAULT_VOLUME_AIR = 0.6;
const double DEFAULT_EPS2 = 0.004 * 0.004;

struct fluid_config {
	double mass_particle, volume_air, volume_liquid, eps2; 

	fluid_config(double mass_particle = DEFAULT_MASS_PARTICLE, double volume_air = DEFAULT_VOLUME_AIR,
				 double eps2 = DEFAULT_EPS2) {
		if (volume_air < 0 || volume_air > 1) {
			throw std::runtime_error("Invalid volume air");
		}
		
		this->eps2 = eps2;
		this->volume_air = volume_air;
		this->volume_liquid = 1 - volume_air;
		this->mass_particle = mass_particle;
	}
};

class Fluid {
private:
    std::vector<Vector2> particles;
    std::vector<Vector2> velocities;
	fluid_config config;

	FluidOptimalTransport *fluid_step(double dt) {
		double default_lambda = 1. / particles.size() * config.volume_liquid;
		std::vector<double> lambdas(particles.size(), default_lambda);
		FluidOptimalTransport *ot_solver = new FluidOptimalTransport(particles, std::move(lambdas), config.volume_air);
		const auto& cells = ot_solver->get_cells();

		for(size_t i = 0; i < particles.size(); ++i) {
			Vector2 gravity = Vector2(0., -G * config.mass_particle);
			Vector2 centroid = cells[i].get_center();
			Vector2 otForce = 1. / config.eps2 * (centroid - particles[i]); 
 			Vector2 F_total = gravity + otForce;
			velocities[i] += dt / config.mass_particle * F_total;
			particles[i] += dt * velocities[i];
		}

		return ot_solver;
	}
public:
	Fluid(size_t N_particles, fluid_config conf = fluid_config()): config(conf) {
		particles.reserve(N_particles);
		velocities.reserve(N_particles);
		for(size_t i = 0; i < N_particles; ++i) {
			particles.emplace_back((rand() / double(RAND_MAX), rand() / double(RAND_MAX)));
			velocities.emplace_back(Vector2(0., 0.));
		}
	}

	template<typename SaveFun>
	void run_sim(size_t steps, double total_time, SaveFun save_frame) {
		double dt = total_time / double(steps);
		for(size_t i = 0; i < steps; ++i) {
			FluidOptimalTransport *ot_solver = fluid_step(dt);
			save_frame(ot_solver->get_cells(), i);
			delete ot_solver;
		}
	}
};