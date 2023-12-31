#pragma once 

#include "vector.h"
#include "poly.hpp"
#include "../benchmarker.hpp"

#include <vector>

#define VORONOI_SIMPLE 0
#define VORONOI_PARALLEL 1

#ifndef VORONOI_ALGORITHM
#define VORONOI_ALGORITHM VORONOI_PARALLEL
#endif

#if VORONOI_ALGORITHM == VORONOI_SIMPLE
class PowerDiagram {
protected:
    using WeightedPoint = std::pair<Vector2, double>;

    std::vector<Vector2> points;
    std::vector<double> weights;

    std::vector<Polygon> cells;

    std::vector<Vector2> clip_by_bisector(const std::vector<Vector2>& vertices, const WeightedPoint& pair_0, const WeightedPoint& pair_i) {
        std::vector<Vector2> result;

        const Vector2& p0 = pair_0.first;
        const Vector2& pi = pair_i.first;
        double w0 = pair_0.second;
        double wi = pair_i.second;

		Vector2 M = (p0 + pi) * 0.5;
		double d2 = (p0 - pi).norm2();	
		Vector2 Mprime = M + (w0 - wi) / (2 * d2) * (pi - p0);
		for(int i = 0; i < vertices.size(); ++i) {
			const Vector2& A = i == 0 ? vertices.back() : vertices[i - 1];
			const Vector2& B = vertices[i];

			// Check if A-B is inside clipping domain or in-between
			double t = (Mprime - A).dot(pi - p0) / (B - A).dot(pi - p0);
			Vector2 P = A + t * (B - A);

            auto dist_0 = [w0, p0](const Vector2& v) {
                return (v - p0).norm2() - w0;
            };

            auto dist_i = [wi, pi](const Vector2& v) {
                return (v - pi).norm2() - wi;
            };

			if (dist_0(B) < dist_i(B)) {
				// B is inside
				if (dist_0(A) >= dist_i(A)) {
					// A is outside
					result.emplace_back(P);
				} 
				result.emplace_back(B);
			}
			else {
				// B is outside
				if (dist_0(A) < dist_i(A)) {
					// A is inside
					result.emplace_back(P);
				}
			}
		}

		return std::move(result);
    }

    virtual void compute_cells() {
        Benchmarker::start_one("compute_cells");
        cells.reserve(points.size());
        // NOTE: I think this makes it slower 
        // #pragma omp parallel for schedule(dynamic) shared(cells)
        for(size_t i = 0; i < points.size(); ++i) {
            std::vector<Vector2> vertices = {  
                Vector2(0, 0),
                Vector2(0, 1),
                Vector2(1, 1),
                Vector2(1, 0)
            };

            for(size_t j = 0; j < points.size(); ++j) {
                if(i == j) 
                    continue;

                vertices = clip_by_bisector(vertices, {points[i], weights[i]}, {points[j], weights[j]});
	    	}

            cells.emplace_back(Polygon(std::move(vertices)));
		}

        Benchmarker::end_one("compute_cells");
    }
public:
    PowerDiagram() = default;

    PowerDiagram(const std::vector<Vector2>& points, const std::vector<double>& weights) {
        if (points.size() != weights.size()) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = points;
        this->weights = weights;
    }

    PowerDiagram(const std::vector<Vector2>& points, double *weights, size_t n) {
        if (points.size() != n) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = points;
        this->weights.resize(n);
        for(size_t i = 0; i < n; ++i) {
            this->weights[i] = weights[i];
        }
    }

    PowerDiagram(std::vector<Vector2>&& points, std::vector<double>&& weights) {
        if (points.size() != weights.size()) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = std::move(points);
        this->weights = std::move(weights);
    }

    void update_weights(const double *x, size_t n) {
        if (n != weights.size()) {
            throw std::runtime_error("Number of weights must be equal to number of points");
        }

        for(size_t i = 0; i < n; ++i) {
            weights[i] = x[i];
        }

        // Invalidate cells so that they are recomputed on next call to get_cells()
        // This is necessary because the cells depend on the weights
        cells.clear();
    }

    const std::vector<Polygon>& get_cells() {
        if (cells.empty())
            compute_cells();

        return cells;
    }

    const std::vector<Vector2>& get_points() const {
        return points;
    }

    const std::vector<double>& get_weights() const {
        return weights;
    }
};
#else 
#include "nanoflann.hpp"
#include "vov_adaptor.h"

class PowerDiagram {
protected:
    using WeightedPoint = std::pair<Vector2, double>;

    std::vector<Vector2> points;
    std::vector<double> weights;

    std::vector<Polygon> cells;

    std::vector<Vector2> clip_by_bisector(const std::vector<Vector2>& vertices, const WeightedPoint& pair_0, const WeightedPoint& pair_i) {
        std::vector<Vector2> result;

        const Vector2& p0 = pair_0.first;
        const Vector2& pi = pair_i.first;
        double w0 = pair_0.second;
        double wi = pair_i.second;

		Vector2 M = (p0 + pi) * 0.5;
		double d2 = (p0 - pi).norm2();	
		Vector2 Mprime = M + (w0 - wi) / (2 * d2) * (pi - p0);
		for(int i = 0; i < vertices.size(); ++i) {
			const Vector2& A = i == 0 ? vertices.back() : vertices[i - 1];
			const Vector2& B = vertices[i];

			// Check if A-B is inside clipping domain or in-between
			double t = (Mprime - A).dot(pi - p0) / (B - A).dot(pi - p0);
			Vector2 P = A + t * (B - A);

            auto dist_0 = [w0, p0](const Vector2& v) {
                return (v - p0).norm2() - w0;
            };

            auto dist_i = [wi, pi](const Vector2& v) {
                return (v - pi).norm2() - wi;
            };

			if (dist_0(B) < dist_i(B)) {
				// B is inside
				if (dist_0(A) >= dist_i(A)) {
					// A is outside
					result.emplace_back(P);
				} 
				result.emplace_back(B);
			}
			else {
				// B is outside
				if (dist_0(A) < dist_i(A)) {
					// A is inside
					result.emplace_back(P);
				}
			}
		}

		return std::move(result);
    }

    virtual void compute_cells() {
        using namespace nanoflann;
        using my_vector_of_vectors_t = std::vector<std::array<double, 3>>;
        using my_kd_tree_t = KDTreeVectorOfVectorsAdaptor<my_vector_of_vectors_t, double>;
        cells.reserve(points.size());
        double max_weight = *std::max_element(weights.begin(), weights.end());
        std::vector<std::array<double, 3>> weighted_points(points.size());
        for(size_t i = 0; i < points.size(); ++i) {
            // Take the std::abs just to make sure we don't take sqrt of -eps
            double last_corrd = std::sqrt(std::abs(max_weight - weights[i]));
            weighted_points[i] = {points[i].x, points[i].y, last_corrd};
        }
        my_kd_tree_t tree(3 /*dim*/, weighted_points, 10 /* max leaf */);
        for(size_t i = 0; i < points.size(); ++i) {
            std::vector<Vector2> vertices = {  
                Vector2(0, 0),
                Vector2(0, 1),
                Vector2(1, 1),
                Vector2(1, 0)
            };

            const size_t k = 10;
            WeightedPoint pair_0 = {points[i], weights[i]};
            double R = 0;
            double max_dist = 0;
            size_t idx = 1;
            do {
                size_t num_results = idx * k + 1;
                // query k nearest neighbors
                std::vector<size_t> ret_indexes(num_results);
                std::vector<double> out_dists_sqr(num_results);
                nanoflann::KNNResultSet<double> resultSet(num_results);
                resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                tree.index->findNeighbors(resultSet, &weighted_points[i][0]);
                for(size_t j = (idx - 1) * k + 1; j < resultSet.size(); ++j) {
                    size_t pts_index = ret_indexes[j];
                    vertices = clip_by_bisector(vertices, pair_0, {points[pts_index], weights[pts_index]});
                }
                max_dist = std::sqrt(out_dists_sqr.back()); // - weights[i];
                R = 0;
                for(auto& vertex : vertices) {
                    double dist2 = points[i].dist2(vertex) + max_weight - weights[i];
                    R = std::max(R, std::sqrt(dist2));
                }
                ++idx;
            } while (2 * R >= max_dist);

            cells.emplace_back(Polygon(std::move(vertices)));
		}
    }
public:
    PowerDiagram() = default;

    PowerDiagram(const std::vector<Vector2>& points, const std::vector<double>& weights) {
        if (points.size() != weights.size()) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = points;
        this->weights = weights;
    }

    PowerDiagram(const std::vector<Vector2>& points, double *weights, size_t n) {
        if (points.size() != n) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = points;
        this->weights.resize(n);
        for(size_t i = 0; i < n; ++i) {
            this->weights[i] = weights[i];
        }
    }

    PowerDiagram(std::vector<Vector2>&& points, std::vector<double>&& weights) {
        if (points.size() != weights.size()) {
            throw std::runtime_error("Number of points and weights must be equal");
        }

        this->points = std::move(points);
        this->weights = std::move(weights);
    }

    void update_weights(const double *x, size_t n) {
        if (n != weights.size()) {
            throw std::runtime_error("Number of weights must be equal to number of points");
        }

        for(size_t i = 0; i < n; ++i) {
            weights[i] = x[i];
        }

        // Invalidate cells so that they are recomputed on next call to get_cells()
        // This is necessary because the cells depend on the weights
        cells.clear();
    }

    const std::vector<Polygon>& get_cells() {
        if (cells.empty())
            compute_cells();

        return cells;
    }

    const std::vector<Vector2>& get_points() const {
        return points;
    }

    const std::vector<double>& get_weights() const {
        return weights;
    }
};
#endif 
