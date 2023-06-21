#pragma once 
#include <vector>
#include <string.h>

#include "lbfgs.h"
#include "vector.h"
#include "diagram.hpp"
#include "poly.hpp"

class FluidPowerDiagram : public PowerDiagram {
private:
    double w_air;
    Polygon disk;

    void compute_unit_disk() {
        const size_t DISK_SIZE = 30;
        // Avoid recomputing the disk vertices every time
        std::vector<Vector2> disk_vertices;
        disk_vertices.reserve(DISK_SIZE);
		for(size_t i = 0; i < DISK_SIZE; ++i) {
			double t = ((double)i / double(DISK_SIZE)) * 2 * M_PI;
			disk_vertices.emplace_back(Vector2(std::cos(t), std::sin(t)));
		}

        disk = Polygon(std::move(disk_vertices));
    }

    void compute_cells() override {
        Benchmarker::start_one("Computing cells");
        PowerDiagram::compute_cells();
        Benchmarker::end_one("Computing cells");
        // Avoid recomputing the cells every time
        Benchmarker::start_one("Clipping cells");
        #pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < points.size(); ++i) {
            double radius = std::sqrt(weights[i] - w_air);
            cells[i] = cells[i].clip_by(disk.shift_and_scale(points[i], radius));
        }
        Benchmarker::end_one("Clipping cells");
    }
public:
    FluidPowerDiagram() = default;

    FluidPowerDiagram(const std::vector<Vector2>& points, const std::vector<double>& weights, double w_air) : PowerDiagram(points, weights), w_air(w_air) {
        compute_unit_disk();
    }

    FluidPowerDiagram(const std::vector<Vector2>& points, double *weights, size_t n, double w_air) : PowerDiagram(points, weights, n), w_air(w_air) {
        compute_unit_disk();
    }

    double get_w_air() const {
        return w_air;
    }

    void set_w_air(double w_air) {
        // Invalidate the cells cache
        this->cells.clear();
        this->w_air = w_air;
    }
};

class FluidOptimalTransport {
protected:
    bool computed = false;
    std::vector<Vector2> points;
    std::vector<double> lambdas;
    double target_volume_air;

    FluidPowerDiagram diagram_solver;

    void solve() {
        size_t n = points.size();
        double *x = new double[n + 1];
        // Sad: memset does not work with doubles
        for(size_t i = 0; i < n; ++i) 
            x[i] = 1.;
        
        // Last parameter: weight of air 
        const double initial_w_air = 0.;
        x[n] = initial_w_air;
        diagram_solver = FluidPowerDiagram(points, x, n, initial_w_air);

		double fx = 0;

        Benchmarker::start_one("L-BFGS");
		int ret = lbfgs(n + 1, x, &fx, _evaluate, _progress, this, NULL);
        Benchmarker::end_one("L-BFGS");
        #ifdef LBFGS_LOGGING
        std::cout << "L-BFGS optimization terminated with status code = " << ret << std::endl;
        if (ret < 0) {
            std::cout << "Error description: " << lbfgs_strerror(ret) << std::endl;
        }
        #endif
        delete x;
        computed = true;
	}

    static lbfgsfloatval_t _evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
    {
        return reinterpret_cast<FluidOptimalTransport*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    )
    {
        lbfgsfloatval_t fx = 0.0;

        size_t particles_cnt = n - 1;
        size_t w_air_index = n - 1;
        double w_air = x[w_air_index];

        diagram_solver.update_weights(x, particles_cnt);
        diagram_solver.set_w_air(w_air);
        const auto& cells = diagram_solver.get_cells();
        const auto& weights = diagram_solver.get_weights();

        Benchmarker::start_one("L-BFGS evaluation");
        double current_volume_fluid = 0.;
        for (size_t i = 0; i < particles_cnt; ++i) {
            // if (weights[i] - w_air < -0.0001) {
            //     std::stringstream ss;
            //     ss << "Particle weight is less than air weight: " << weights[i] << " < " << w_air << " at index " << i;
            //     throw std::runtime_error(ss.str());
            // }

			double cell_area = cells[i].area();
            current_volume_fluid += cell_area;
			g[i] = -(lambdas[i] - cell_area);
			// first term, int |x - y|^2 f(x) dx
			fx += cells[i].int_norm_2(points[i]);
			// second term, weight * area 
			fx -= cell_area * weights[i];
			// third term, lambda * weight
			fx += lambdas[i] * weights[i];
        }

        // if (current_volume_fluid > 1.0001) {
        //     throw std::runtime_error("Fluid volume is greater than 1: " + std::to_string(current_volume_fluid));
        // }

        double current_volume_air = 1. - current_volume_fluid;
        double vol_difference = target_volume_air - current_volume_air;
        g[w_air_index] = -vol_difference;
        fx += w_air * vol_difference;

        Benchmarker::end_one("L-BFGS evaluation");
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<FluidOptimalTransport*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        #ifdef LBFGS_LOGGING
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        #endif
        return 0;
    }
public:
    FluidOptimalTransport(const std::vector<Vector2>& points, const std::vector<double>& lambdas, double target_volume_air) : points(points), lambdas(lambdas) {
        if (points.size() != lambdas.size()) {
            throw std::runtime_error("Number of points and lambdas must be equal");
        }

        this->target_volume_air = target_volume_air;
    }

    FluidOptimalTransport(const std::vector<Vector2>& points, std::vector<double>&& lambdas, double target_volume_air) {
        if (points.size() != lambdas.size()) {
            throw std::runtime_error("Number of points and lambdas must be equal");
        }

        this->points = points;
        this->lambdas = std::move(lambdas);

        this->target_volume_air = target_volume_air;
    }

    const std::vector<Polygon>& get_cells() {
        if (!computed)
            solve();

        return diagram_solver.get_cells();
    }
};