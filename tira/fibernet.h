#pragma once

#include "fiber.h"

namespace tira {

    /**
     * @brief      Class provides functionality for an interconnected network of fibers.
     *
     * @tparam     VertexAttribute  The data type used to store user-defined attributes at each vertex
     * @tparam     EdgeAttribute    The data type used to store user-defined attributes at each edge
     * @tparam     NodeAttribute    The data type used to store user-defined attributes at each node
     */
	template <typename VertexAttribute = float, 
        typename EdgeAttribute = unsigned int,
        typename NodeAttribute = unsigned int>
	class fibernet {

    public:

        /**
         * @brief      This edge class extends a fiber to include a user-defined edge attribute 
         * and a list of connected nodes
         */
        class edge : public fiber<VertexAttribute> {
        protected:
            EdgeAttribute _ea;
            size_t _n[2];              // node indices

        public:

            /**
             * @brief      Constructs a new edge given an existing fiber, two nodes, and an edge attribute
             *
             * @param[in]  f       existing fiber containing a series of geometric vertices and their attributes
             * @param[in]  n0      first node connected by this edge (closest to the first vertex in the fiber)
             * @param[in]  n1      second node connected to this edge (closest to the last vertex in the fiber)
             * @param[in]  attrib  user-defined attribute to store information within each edge
             */
            edge(fiber<VertexAttribute> f, size_t n0, size_t n1, EdgeAttribute attrib) : fiber<VertexAttribute>(f) {
                _n[0] = n0;
                _n[1] = n1;
                _ea = attrib;
            }

            /**
             * @brief      Creates an empty edge linking two nodes
             *
             * @param[in]  n0      first node ID connected by this edge
             * @param[in]  n1      second node ID connected by this edge
             * @param[in]  attrib  user-defined attribute that stores information at this edge
             */
            edge(size_t n0, size_t n1, EdgeAttribute attrib) {
                _n[0] = n0;
                _n[1] = n1;
                _ea = attrib;
            }

            /**
             * @brief      Assigns or adjusts nodes connected by this edge.
             *
             * @param[in]  n0    sets the first node connected by this edge
             * @param[in]  n1    sets the second node connected by this edge
             */
            void nodes(size_t n0, size_t n1) {
                _n[0] = n0;
                _n[1] = n1;
            }

            /**
             * @brief      Retrieve the first node ID
             *
             * @return     index into the _nodes vector
             */
            size_t n0() { return _n[0]; }

            /**
             * @brief      Retrieve the first node ID
             *
             * @return     index into the _nodes vector
             */
            size_t n1() { return _n[1]; }

            /**
             * @brief      Returns the user-specified edge attribute
             *
             * @return     The edge attribute.
             */
            EdgeAttribute ea() { return _ea; }

            /**
             * @brief      Assigns an attribute to this edge
             *
             * @param[in]  attrib  The attribute
             */
            void ea(EdgeAttribute attrib) { _ea = attrib; }
          

        };

        /**
         * @brief      Defines a class for a node, which connects multiple fibers and consists of a geometric point and attributes
         */
        class node : public vertex<VertexAttribute> {
        protected:
            NodeAttribute _na;
            std::vector<size_t> _ei;   // indices of connected edges

        public:

            /**
             * @brief      Create a new node from an existing vertex and a node-specific attribute
             *
             * @param[in]  v       Existing vertex
             * @param[in]  attrib  Node attribute
             */
            node(vertex<VertexAttribute> v, NodeAttribute attrib) : vertex<VertexAttribute>(v) {
                _na = attrib;
            }

            /**
             * @brief      Retrieves the indices of all edges connected to this node
             *
             * @return     vector of index values corresponding to each connected edge
             */
            std::vector<size_t> edges() { return _ei; }

            /**
             * @brief      Connects this node to an edge. This should be done in parity with connecting edges to nodes
             *
             * @param[in]  ei    ID of the edge connecting to this node
             */
            void add_edge(size_t ei) { _ei.push_back(ei); }

            /**
             * @brief      Retrieve the node-specific attribute
             *
             * @return     The node attribute.
             */
            NodeAttribute na() { return _na; }

            /**
             * @brief      Assign or update an attribute for this node
             *
             * @param[in]  na    The new value for the node attribute
             */
            void na(NodeAttribute na) { _na = na; }
        };

    protected:
        /**
         * List of nodes in this network, essentially modeled as a graph
         */
        std::vector<node> _nodes;

        /**
         * List of edges in this network, modeled as a mathematical graph
         */
        std::vector<edge> _edges;

        /**
         * An axis-aligned bounding box containing all geometric positions in the network.
         */
        std::pair<glm::vec3, glm::vec3> _aabb; // bounding box around all fiber points

    public:

        /**
         * @brief      Adds a node to the network using a geometric position and attribute
         *
         * @param[in]  v     vertex describing the geometric position of the node
         * @param[in]  n     node-specific attributes
         *
         * @return     { description_of_the_return_value }
         */
        size_t add_node(vertex<VertexAttribute> v, NodeAttribute n) {
            node new_node(v, n);
            _nodes.push_back(new_node);
            return _nodes.size() - 1;
        }

        /**
         * Add an edge to the graph given its connected nodes and fiber geometry
         * @param inode0 the first node in the graph (associated with the medial axis point nearest the first point in pts)
         * @param inode1 the second node in the graph (associated with the medial axis point closest to the last point in pts)
         * @param f is the fiber geometry
         * @return the internal ID of the edge in the graph
         */
        size_t add_edge(size_t inode0, size_t inode1, fiber<VertexAttribute> f, EdgeAttribute a) {

            edge new_edge(f, inode0, inode1, a);                                  // create a new edge structure
            _edges.push_back(new_edge);
            
            size_t idx = _edges.size() - 1;
            _nodes[inode0].add_edge(idx);
            _nodes[inode1].add_edge(idx);

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

            v0 = _nodes[n0].v;
            v1 = _nodes[n1].v;
        }

        inline fiber<VertexAttribute> fiber_edge(const size_t id) const {
            return _edges[id];
        }
	};
}
