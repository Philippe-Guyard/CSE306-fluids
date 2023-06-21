#include <iostream>
#include <map>
#include <vector>

#include "vector.h"
#include "mesh_reader.hpp"

class Edge {
public:
    int a, b;
    Edge(int u = -1, int v = -1): a(u), b(v) {}
};

bool operator<(const Edge& e1, const Edge& e2) {
    if (std::min(e1.a, e1.b) < std::min(e2.a, e2.b))
        return true; 
    
    return std::max(e1.a, e1.b) < std::max(e1.a, e2.b);
}


int main() {
    // Read goethe.obj with the mesh reader
    TriangleMeshDescriptor mesh;
    mesh.readOBJ("goethe.obj");

    std::map<Edge, std::vector<size_t>> edge_to_tri;
    std::map<int, std::vector<int>> vertex_to_tri;
    for(int i = 0; i < mesh.indices.size(); ++i) {
        int a = mesh.indices[i].vtxi;
        int b = mesh.indices[i].vtxj;
        int c = mesh.indices[i].vtxk;

        edge_to_tri[Edge(a, b)].push_back(i);
        edge_to_tri[Edge(b, c)].push_back(i);
        edge_to_tri[Edge(c, a)].push_back(i);

        vertex_to_tri[a].push_back(i);
        vertex_to_tri[b].push_back(i);
        vertex_to_tri[c].push_back(i);
    }

    std::cout << "Lol" << std::endl;

    std::vector<Edge> boundary_edges;
    std::vector<bool> is_boundary(mesh.vertices.size(), false);
    for(auto it = edge_to_tri.begin(); it != edge_to_tri.end(); it++) {
        if (it->second.size() == 1) {
            boundary_edges.push_back(it->first);
            is_boundary[it->first.a] = true;
            is_boundary[it->first.b] = true;
        }
    }

    std::cout << "Lol 2" << std::endl;

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    ordered_boundary_edges[0] = boundary_edges[0];
    for(int i = 1; i < boundary_edges.size(); ++i) {
        for(int j = 0; j < boundary_edges.size(); ++j) {
            if (boundary_edges[j].a == ordered_boundary_edges[i - 1].b) {
                ordered_boundary_edges[i] = boundary_edges[j];
                break;
            }
        }
    }

    std::cout << "Lol 3" << std::endl;

    for(int i = 0; i < ordered_boundary_edges.size(); ++i) {
        double theta = i / (double)ordered_boundary_edges.size() * 2 * M_PI;
        Vector3 circle_vtx(0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta), 0.);        
        mesh.vertices[ordered_boundary_edges[i].a] = circle_vtx;
    }

    std::cout << "Lol 4" << std::endl;

    for(int iter = 0; iter < 10; ++iter) {
        std::vector<Vector3> updated_vertices = mesh.vertices;
        for(int i = 0; i < mesh.vertices.size(); ++i) {
            if (is_boundary[i])
                continue;

            Vector3 sum_neighbors(0., 0., 0);
            size_t total_neighbors = 0;
            for(int j = 0; j < vertex_to_tri[i].size(); ++j) {
                size_t tri_idx = vertex_to_tri[i][j];
                if (mesh.indices[tri_idx].vtxi != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxi];
                if (mesh.indices[tri_idx].vtxj != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxj];
                if (mesh.indices[tri_idx].vtxk != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxk];

                total_neighbors += 2;
            }
            updated_vertices[i] = sum_neighbors / total_neighbors;
        }

        mesh.vertices = updated_vertices;
    }

    std::cout << "Lol 5" << std::endl;

    //mesh.writeOBJ("goethe_new.obj");

    return 0;
}