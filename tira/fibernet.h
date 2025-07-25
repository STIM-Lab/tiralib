#pragma once

#include "fiber.h"

namespace tira {

	template <typename VertexAttribute = float, 
        typename EdgeAttribute = unsigned int,
        typename NodeAttribute = unsigned int>
	class fibernet {

    protected:
        struct _Edge : public fiber<VertexAttribute> {
            EdgeAttribute a;
            size_t inodes[2];              // node indices

            _Edge(fiber<VertexAttribute> f, size_t n0, size_t n1, EdgeAttribute attrib) : fiber<VertexAttribute>(f) {
                inodes[0] = n0;
                inodes[1] = n1;
                a = attrib;
            }
        };

        struct _Node : public vertex<VertexAttribute> {
            NodeAttribute a;
            std::vector<size_t> iedges;   // indices of connected edges

            _Node(vertex<VertexAttribute> v, NodeAttribute attrib) : vertex<VertexAttribute>(v) {
                a = attrib;
            }
        };

        std::vector<_Node> _nodes;
        std::vector<_Edge> _edges;

        std::pair<glm::vec3, glm::vec3> _aabb; // bounding box around all fiber points

    public:

        size_t add_node(vertex<VertexAttribute> v, NodeAttribute n) {
            _Node new_node(v, n);
            _nodes.push_back(new_node);
        }

        /**
         * Add an edge to the graph given its connected nodes and fiber geometry
         * @param inode0 the first node in the graph (associated with the medial axis point nearest the first point in pts)
         * @param inode1 the second node in the graph (associated with the medial axis point closest to the last point in pts)
         * @param f is the fiber geometry
         * @return the internal ID of the edge in the graph
         */
        size_t add_edge(size_t inode0, size_t inode1, fiber<VertexAttribute> f, EdgeAttribute a) {

            _Edge new_edge(f, inode0, inode1, a);                                  // create a new edge structure
            //new_edge = f;
            //new_edge.inodes[0] = inode0;                    // set the two node IDs based on user input
            //new_edge.inodes[1] = inode1;
            
            //new_edge.ipts.resize(pts.size());               // resize the medial axis point index array to match the number of points in the edge
            //std::iota(new_edge.ipts.begin(), new_edge.ipts.end(), _points.size());  // create a sequential list of IDs starting with the last index in _points
            _edges.push_back(new_edge);

            //_points.insert(_points.end(), pts.begin(), pts.end());      // insert all of the new points for this edge into the _points array

            return _edges.size() - 1;
        }

        /**
     * Returns the 3D coordinates and radii associated with both nodes associated with an edge.
     * @param id identifier for the edge that will be returned
     * @return std::pair storing the coordinates and radii of both nodes connected by the edge
     */
        void graph_edge(const size_t id, vertex<VertexAttribute>& v0, vertex<VertexAttribute>& v1) const {
            const size_t n0 = _edges[id].inodes[0];
            const size_t n1 = _edges[id].inodes[1];

            //std::pair<glm::vec4, glm::vec4> edge_pts;
            //edge_pts.first = _points[n0];
            //edge_pts.second = _points[n1];
            v0 = _nodes[n0].v;
            v1 = _nodes[n1].v;
        }

        inline fiber<VertexAttribute> fiber_edge(const size_t id) const {
            return _edges[id];
        }
	};
}
