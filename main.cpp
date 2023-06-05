#include <vector>
#include <fstream>
#include <iostream>
#include <optional>
#include "lbfgs.h"
#include <chrono>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include "stb_image.h"
#include <sstream>

#include "fluids/vector.h"
#include "fluids/poly.hpp"
#include "fluids/diagram.hpp"

using Vector = Vector2;

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none", std::optional<std::vector<Vector2>> points = std::nullopt) {
	FILE* f = fopen(filename.c_str(), "w+"); 
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (int i=0; i<polygons.size(); i++) {
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \""); 
		for (int j = 0; j < polygons[i].get_size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertex_at(j).x * 1000), (1000 - polygons[i].vertex_at(j).y * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
	}
	if (points.has_value()) {
		// Draw each point as a circle with radius 5
		for (int i = 0; i < points.value().size(); ++i) {
			fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"5\" stroke = \"black\" fill = \"black\"/>\n", points.value()[i].x * 1000, (1 - points.value()[i].y) * 1000);
		}
	}
	fprintf(f, "</svg>\n");
	fclose(f);
}

int sgn(double x) {
	if (x > 0)
		return 1.;
	if (x < 0)
		return -1;
	
	return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
	int W = 500, H = 500;
	std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < cells.size(); i++) {

		double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
		for (int j = 0; j < cells[i].get_size(); j++) {
			bminx = std::min(bminx, cells[i].vertex_at(j).x);
			bminy = std::min(bminy, cells[i].vertex_at(j).y);
			bmaxx = std::max(bmaxx, cells[i].vertex_at(j).x);
			bmaxy = std::max(bmaxy, cells[i].vertex_at(j).y);
		}
		bminx = std::min(W-1., std::max(0., W * bminx));
		bminy = std::min(H-1., std::max(0., H * bminy));
		bmaxx = std::max(W-1., std::max(0., W * bmaxx));
		bmaxy = std::max(H-1., std::max(0., H * bmaxy));

		for (int y = bminy; y < bmaxy; y++) {
			for (int x = bminx; x < bmaxx; x++) {
				int prevSign = 0;
				bool isInside = true;
				double mindistEdge = 1E9;
				for (int j = 0; j < cells[i].get_size(); j++) {
					double x0 = cells[i].vertex_at(j).x * W;
					double y0 = cells[i].vertex_at(j).y * H;
					double x1 = cells[i].vertex_at((j + 1) % cells[i].get_size()).x * W;
					double y1 = cells[i].vertex_at((j + 1) % cells[i].get_size()).y * H;
					double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
					int sign = sgn(det);
					if (prevSign == 0) prevSign = sign; else
						if (sign == 0) sign = prevSign; else
						if (sign != prevSign) {
							isInside = false;
							break;
						}
					prevSign = sign;
					double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
					double distEdge = std::abs(det)/ edgeLen;
					double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
					if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
					mindistEdge = std::min(mindistEdge, distEdge);
				}
				if (isInside) {
					//if (i < N) {   // the N first particles may represent fluid, displayed in blue
					//	image[((H - y - 1)*W + x) * 3] = 0;
					//	image[((H - y - 1)*W + x) * 3 + 1] = 0;
					//	image[((H - y - 1)*W + x) * 3 + 2] = 255;
					//}
					if (mindistEdge <= 2) {
						image[((H - y - 1)*W + x) * 3] = 0;
						image[((H - y - 1)*W + x) * 3 + 1] = 0;
						image[((H - y - 1)*W + x) * 3 + 2] = 0;
					}

				}
				
			}
		}
	}
	std::ostringstream os;
	os << filename << frameid << ".png";
	stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

void save(PowerDiagram& pd, std::string filename) {
	std::vector<Polygon> cells = pd.get_cells();
	std::vector<Vector2> points = pd.get_points();
	save_svg(cells, filename, "none", points);
}

int main() {
	srand(time(0));
	std::vector<Vector> points = {
        Vector(0.23, 0.54),
        Vector(0.63, 0.23),
        Vector(0.34, 0.31),
        Vector(0.45, 0.67),
        Vector(0.56, 0.45),
        Vector(0.67, 0.56),
        Vector(0.78, 0.67),
        Vector(0.89, 0.78),
        Vector(0.12, 0.89),
    };

	std::vector<double> weights;
	for (int i = 0; i < points.size(); ++i) {
		weights.emplace_back(1. / (double)points.size());
	}

	PowerDiagram pd(points, weights);
	save(pd, "voronoi.svg");

	return 0;
}