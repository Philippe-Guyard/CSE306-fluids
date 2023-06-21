#pragma once 

#include <vector>
#include <algorithm>
#include <queue>

#include "vector.h"

size_t k = 2;
using TVec = Vector2;

// TODO: This should be templated, but couldn't be bothered 
// template <typename TVec, std::size_t k>
class KdTree {
using iter = std::vector<TVec>::iterator;

private: 
    struct node {
        TVec location, bbox_min, bbox_max;
        node *left, *right; 

        node(const TVec& loc, const TVec& low, const TVec& high): location(loc), bbox_min(low), bbox_max(high) {}

        ~node() {
            if (left != nullptr)
                delete left;
            if (right != nullptr)
                delete right; 
        }

        // Check if box contains point
        bool bbox_contains(const TVec& p) {
            for(size_t j = 0; j < k; ++j) {
                if (p[j] < bbox_min[j] || p[j] > bbox_max[j])
                    return false;
            }

            return true;
        }

        // Check if box contains sphere 
        bool bbox_contains(const TVec& p, double r) {
            for(size_t j = 0; j < k; ++j) {
                if ((p[j] - bbox_min[j]) * (p[j] - bbox_min[j]) > r * r)
                    return false;
            }

            return true;
        }
    };

    size_t size;
    KdTree::node *root;

    KdTree::node* build(iter begin, iter end, int depth,
                        const TVec& low, const TVec& high) {
        if (end == begin)
            return nullptr;

        if (end == begin + 1) {
            return new node(*begin, low, high);
        }
        size_t m = (end - begin) / 2;
        int axis = depth % k;
        std::nth_element(begin, begin + m, end, [axis](const TVec& a, const TVec& b) {
            return a[axis] < b[axis];
        });

        auto node_loc = *(begin + m);
        KdTree::node* result = new node(node_loc, low, high);
        // No need to waste time recomputing a strict bounding box for children.
        // In the worst case (highly probably), only the values along "axis"
        // will change, so we just update that 
        if (m > 0) {
            TVec new_high = high;
            new_high[axis] = node_loc[axis];
            result->left = build(begin, begin + m, depth + 1, low, new_high);
        }
        if (begin + m + 1 != end) {
            TVec new_low = low;
            new_low[axis] = node_loc[axis];
            result->right = build(begin + m + 1, end, depth + 1, new_low, high);
        }
    }
public:
    // Copy the argument 
    KdTree(std::vector<TVec> points) {
        TVec points_min = points[0], points_max = points[0];
        for(size_t i = 1; i < points.size(); ++i) {
            for(size_t j = 0; j < k; ++j) {
                points_max[j] = std::max(points_max[j], points[i][j]);
                points_min[j] = std::min(points_min[j], points[i][j]);
            }
        }

        this->root = build(points.begin(), points.end(), 0, points_min, points_max);
        this->size = points.size();
    }

    ~KdTree() {
        delete root;
    }

    std::vector<TVec> n_nearest(const TVec& source) {
        using temp = std::pair<double, node*>;
        std::priority_queue<temp>();
    }
};