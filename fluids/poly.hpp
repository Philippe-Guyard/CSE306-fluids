#pragma once 

#include <iostream>
#include <vector>
#include <optional>
#include <array>

#include "vector.h"

class Polygon {  
private:
    std::vector<Vector2> vertices;
public:
    Polygon(std::vector<Vector2>&& v) : vertices(std::move(v)) { }

    Polygon(const std::vector<Vector2>& v) : vertices(v) { }

    Polygon() = default;

	double area() const {
		if (vertices.size() < 3)
			return 0;

		double result = 0;
		for (int i = 0; i < vertices.size(); ++i) {
			const Vector2& v1 = vertices[i];
			const Vector2& v2 = vertices[(i + 1) % vertices.size()];
			result += v1.x * v2.y - v1.y * v2.x;
		}

		#ifdef WARNINGS
		if (std::fabs(result) < 1e-7) {
			std::cout << "Warning: area of polygon is too small: " << result << std::endl;
		}
		#endif 
		return std::fabs(result / 2.);
	}

    const Vector2& vertex_at(size_t i) const {
		const auto& N = vertices.size();
        if (i >= N)
            i -= N;
        if (i < 0)
            i += N;
        
        return vertices[i];
    }

    size_t get_size() const {
        return vertices.size();
    }

	double int_norm_2(const Vector2& p) const {
		if (vertices.size() < 3)
			return 0;

		double value = 0.;
		int num_triangles = vertices.size() - 2;
		for(size_t i = 1; i <= num_triangles; ++i) {
			Vector2 triangles[3] = {vertices[0], vertices[i], vertices[i + 1]};
			double local_value = 0.;
			for(size_t k = 0; k < 3; ++k) {
				for(size_t l = k; l < 3; ++l) {
					local_value += (triangles[k] - p).dot(triangles[l] - p);
				}
			}

			Vector2 e1 = triangles[1] - triangles[0];
			Vector2 e2 = triangles[2] - triangles[0];
			double triangle_area = 0.5 * std::fabs(e1.y * e2.x - e1.x * e2.y);
			if (triangle_area < 1e-7) {
                #ifdef WARNINGS
				std::cout << "Warning: triangle area is too small: " << triangle_area << std::endl;
				std::cout << "e1: " << e1.x << " " << e1.y << std::endl;
				std::cout << "e2: " << e2.x << " " << e2.y << std::endl;
				std::cout << e1.y * e2.x - e1.x * e2.y << std::endl;
                #endif 
				continue;
			}
			value += local_value / 6. * triangle_area;
		}

		return value;
	}

	Vector2 get_center() const {
		double signed_area = 0.;
		Vector2 c(0., 0.);
		size_t N = vertices.size();
		for(int i = 0; i < N; ++i) {
			const Vector2& v1 = vertices[i];
			const Vector2& v2 = vertices[(i + 1) % N];
			double cross = v1.x * v2.y - v2.x * v1.y;
			c.x += (v1.x + v2.x) * cross;
			c.y += (v1.y + v2.y) * cross; 
			signed_area += cross;
		}

		return c / (3 * signed_area);
	}

    Polygon shift_and_scale(const Vector2& v, double s) const {
        std::vector<Vector2> new_vertices;
        for(auto& vertex : vertices) {
            new_vertices.emplace_back(v + vertex * s);
        }

        return Polygon(std::move(new_vertices));
    }

    Polygon clip_by_edge(const Vector2& u, const Vector2& v) {
        std::vector<Vector2> result;
		Vector2 N(v.y - u.y, u.x - v.x);
		for(int i = 0; i < vertices.size(); ++i) {
			// const Vector2& A = vertex_at(i - 1);
			// const Vector2& B = vertex_at(i);
			const Vector2& A = i == 0 ? vertices.back() : vertices[i - 1];
			const Vector2& B = vertices[i];
			// Check if A-B is inside clipping domain or in-between
			double t = (u - A).dot(N) / (B - A).dot(N);
			Vector2 P = A + t * (B - A);
			if ((B - u).dot(N) < 0.) {
				// B is inside
				if ((A - u).dot(N) > 0.) {
					// A is outside
					result.emplace_back(P);
				} 
				result.emplace_back(B);
			}
			else {
				// B is outside
				if ((A - u).dot(N) < 0.) {
					// A is inside
					result.emplace_back(P);
				}
			}
		}

		return Polygon(std::move(result));
    }

    Polygon clip_by(const Polygon& other) {
        Polygon intersected = Polygon(vertices); 
		for(size_t i = 0; i < other.vertices.size(); ++i) {
			const Vector2& u = other.vertices[i];
			const Vector2& v = other.vertices[(i + 1) % other.vertices.size()];
			intersected = intersected.clip_by_edge(u, v);
		}

		return intersected;
    }
};	