#pragma once 
#include <vector>
#include <string.h>

#include "lbfgs.h"
#include "vector.h"
#include "diagram.hpp"
#include "poly.hpp"

class OptimalTransport {
private:
    bool computed = false;
    std::vector<Vector2> points;
    std::vector<double> lambdas;

    PowerDiagram diagram_solver;

    void solve() {
        size_t n = points.size();
        double *x = new double[n];
        memset(x, 1., sizeof(double) * n);
        diagram_solver = PowerDiagram(points, x, n);

		double fx = 0;
 
		int ret = lbfgs(n, x, &fx, _evaluate, _progress, this, NULL);
        #ifdef LBFGS_LOGGING
        std::cout << "L-BFGS optimization terminated with status code = " << ret << std::endl;
        #endif
	}

        static lbfgsfloatval_t _evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
    {
        return reinterpret_cast<OptimalTransport*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;

        diagram_solver.update_weights(x, n);
        const auto& cells = diagram_solver.get_cells();
        const auto& weights = diagram_solver.get_weights();

        for (size_t i = 0; i < n; ++i) {
			double cell_area = cells[i].area();
			g[i] = -(lambdas[i] - cell_area);
			// first term, int |x - y|^2 f(x) dx
			fx += cells[i].int_norm_2(points[i]);
			// second term, weight * area 
			fx -= cell_area * weights[i];
			// third term, lambda * weight
			fx += lambdas[i] * weights[i];
        }

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
        return reinterpret_cast<OptimalTransport*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
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
    OptimalTransport(const std::vector<Vector2>& points, const std::vector<double>& lambdas) : points(points), lambdas(lambdas) {
        if (points.size() != lambdas.size()) {
            throw std::runtime_error("Number of points and lambdas must be equal");
        }
    }

    const std::vector<Polygon>& get_cells() {
        if (!computed)
            solve();

        return diagram_solver.get_cells();
    }
};