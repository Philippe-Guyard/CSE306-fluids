#pragma once 

#include <iostream>
#include <vector>
#include <optional>
#include <array>

#include "vector.h"

class Polygon {  
private: 
    size_t N;
	std::vector<Vector2> vertices;

    bool check_vertices() const {
        return N >= 3;
    }

    void assert_vertices() const {
        if (!check_vertices()) {
            #if defined(ERROR)
            throw std::runtime_error("Invalid vertices for polygon");
            #elif defined(WARNINGS)
            std::cout << "Warning: invalid vertices for polygon" << std::endl;
            #endif
        }
    }
public:
    Polygon() = default;

    Polygon(const std::vector<Vector2>& v) : vertices(v), N(v.size()) {
        assert_vertices();
    }

    Polygon(std::vector<Vector2>&& v) {
        vertices = std::move(v);
        N = vertices.size();
        assert_vertices();
    }

    const Vector2& vertex_at(size_t index) const {
        // Support indices from (-N + 1) to (2 * N - 1)
        if (index >= N)
            index -= N;
        if (index < 0)
            index += N;

        return vertices[index];
    }

    size_t get_size() const {
        return N;
    }

    double signed_area() const { 
        if (!check_vertices())
            return 0.;

        double result = 0.;
        for(size_t i = 0; i < N; ++i) {
            const Vector2& v1 = vertex_at(i);
            const Vector2& v2 = vertex_at(i + 1);
            result += v1.cross(v2);
        }
        if (std::fabs(result) < 1e-7) {
            #ifdef WARNINGS
			std::cout << "Warning: area of polygon is too small: " << result << std::endl;
            #endif 
		}

        result *= 0.5;
        return result;
    }

	double area() const {
        return std::fabs(signed_area());
	}

	double int_norm_2(const Vector2& p) const {
        if (!check_vertices())
            return 0.;

        // Compute integral |x - p|^2 f(x) dx over our polygon 
		double result = 0.;
		size_t num_triangles = N - 2;
		for(size_t i = 1; i <= num_triangles; ++i) {
			Vector2 triangles[3] = {vertex_at(0), vertex_at(i), vertex_at(i + 1)};
			double local_value = 0.;
			for(int k = 0; k < 3; ++k) 
				for(int l = k; l < 3; ++l) 
					local_value += (triangles[k] - p).dot(triangles[l] - p);
				
			Vector2 e1 = triangles[1] - triangles[0];
			Vector2 e2 = triangles[2] - triangles[0];
			double triangle_area = 0.5 * std::fabs(e1.cross(e2));
			if (triangle_area < 1e-7) {
                #ifdef WARNINGS
				std::cout << "Warning: triangle area is too small: " << triangle_area << std::endl;
				std::cout << "e1: " << e1.x << " " << e1.y << std::endl;
				std::cout << "e2: " << e2.x << " " << e2.y << std::endl;
				std::cout << e1.y * e2.x - e1.x * e2.y << std::endl;
                #endif
				continue;
			}
			result += local_value / 6. * triangle_area;
		}

		return result;
	}

	Vector2 get_center() const {
        if (!check_vertices())
            return Vector2(0., 0.);
            
		Vector2 c(0., 0.);
		for(int i = 0; i < N; ++i) {
			const Vector2& v1 = vertex_at(i);
			const Vector2& v2 = vertex_at(i + 1);
			double cross = v1.cross(v2);
			c.x += (v1.x + v2.x) * cross;
			c.y += (v1.y + v2.y) * cross; 
		}

		return c / (6 * signed_area());
	}

    Polygon clip_by(const Polygon& other) {
        Polygon result = other;
        for(size_t i = 0; i < other.get_size(); ++i) {
            const Vector2& u = other.vertex_at(i);
            const Vector2& v = other.vertex_at(i + 1);
            Vector2 N(v.y - u.y, u.x - v.x);
            std::vector<Vector2> new_vertices;
            for(size_t j = 0; j < result.get_size(); ++j) {
                const Vector2& A = result.vertex_at(j);
                const Vector2& B = result.vertex_at(j + 1);
                double dA = N.dot(A - u);
                double dB = N.dot(B - u);
                double t = -dA / (B - A).dot(N);
                Vector2 intersection = A + t * (B - A);
                if (dA * dB < 0.) 
                    new_vertices.emplace_back(intersection);
                if (dB < 0.) 
                    new_vertices.emplace_back(B);
            }

            result = Polygon(std::move(new_vertices));
        }

        return result;
    }

    Polygon shift_by(const Vector2& v) const {
        std::vector<Vector2> new_vertices;
        for(size_t i = 0; i < N; ++i) {
            new_vertices.emplace_back(vertex_at(i) + v);
        }

        return Polygon(std::move(new_vertices));
    }

    Polygon scale_by(double s) const {
        std::vector<Vector2> new_vertices;
        for(size_t i = 0; i < N; ++i) {
            new_vertices.emplace_back(vertex_at(i) * s);
        }

        return Polygon(std::move(new_vertices));
    }

    Polygon shift_and_scale(const Vector2& v, double s) const {
        std::vector<Vector2> new_vertices;
        for(size_t i = 0; i < N; ++i) {
            new_vertices.emplace_back((vertex_at(i) + v) * s);
        }

        return Polygon(std::move(new_vertices));
    }
};	