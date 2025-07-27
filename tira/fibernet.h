#pragma once

#include "fiber.h"

namespace tira {

    /**
     * A fibernet object is an undirected spatial graph with the following characteristics:
     * - Spatial locations in the graph are represented by a vertex consisting of a 3D position and
     *    a set of user-specified attributes (indicated by the VertexAttribute template parameter).
     * - Edges are represented by a geometric path between nodes represented by an array of vertex data
     * - Nodes are represented by a single vertex
     * - Edges and Nodes can contain additional independent attributes specified by the EdgeAttribute
     *   and NodeAttribute template parameters
     *
     * IMPORTANT: The graph is designed to be undirected, however the nodes for each edge must be specified
     *   in the correct order so that they align with the path geometry for each edge. Calculations performed
     *   in member functions consider the graph undirected. However, those calculations require the correct
     *   alignment between nodes and path vertices.
     *
     *   For example, consider the following sequence of spatial positions:
     *   [n0] p0---p1---p2---...---pn [n1].
     *   Calculating the length() of this edge requires that n0 be adjacent to p0. In short, Nodes cannot
     *   be trivially swapped within an edge and care should be taken to specify node positions in the
     *   correct order when adding fibers.
     *
     *
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
         * and a list of connected edges
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
            edge(fiber<VertexAttribute> f, size_t n0, size_t n1, EdgeAttribute attrib = {}) : fiber<VertexAttribute>(f) {
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
            edge(size_t n0, size_t n1, EdgeAttribute attrib = {}) {
                _n[0] = n0;
                _n[1] = n1;
                _ea = attrib;
            }

            /**
             * @brief       Creates a new edge using an existing edge and new fiber - Used to update the geometry of the edge.
             *
             * @param f     New fiber geometry
             * @param e     Existing edge, where the graph and attribute parameters will be copied
             */
            edge(fiber<VertexAttribute> f, edge& e) : fiber<VertexAttribute>(f) {
                _n[0] = e._n[0];
                _n[1] = e._n[1];
                _ea = e._ea;
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
            EdgeAttribute& ea() { return _ea; }

            /**
             * @brief      Assigns an attribute to this edge
             *
             * @param[in]  attrib  The attribute
             */
            void ea(EdgeAttribute attrib) { _ea = attrib; }

            /**
             * Smoothing an edge is a little more complicated than smoothing a fiber: we have to account for the node positions,
             * which are not included in the "fiber" component of the edge. We do that by:
             * 1) Create a new fiber that duplicates the current edge fiber
             * 2) Insert both node points at either end of the fiber
             * 3) Smooth the new fiber (while anchoring the end points)
             * 4) Removing the node points from the fiber
             * 5) Creating a new edge from the internal nodes of the smoothed fiber
             *
             * @param sigma the standard deviation of the smoothing kernel (in the position units of the vertices)
             * @return a new edge with a smoothed fiber component
             */
            edge smooth_gaussian(float sigma, vertex<VertexAttribute> node_v0, vertex<VertexAttribute> node_v1) {
                fiber<VertexAttribute> original_fiber = *this;

                original_fiber.insert(original_fiber.begin(), node_v0);
                original_fiber.push_back(node_v1);
                fiber<VertexAttribute> smoothed_fiber = original_fiber.smooth_gaussian(sigma);
                smoothed_fiber.erase(smoothed_fiber.begin());
                smoothed_fiber.pop_back();

                edge smoothed_edge(smoothed_fiber, _n[0], _n[1]);
                return smoothed_edge;
            }

            float length(vertex<VertexAttribute> v0, vertex<VertexAttribute> v1) {
                fiber<VertexAttribute> original_fiber = *this;
                original_fiber.insert(original_fiber.begin(), v0);
                original_fiber.push_back(v1);
                return original_fiber.length();
            }
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

    public:
        /**
         * @brief Returns the vertex representing the "first" node in an edge
         * @param ei ID of the edge
         * @return a vertex<VertexAttribute> providing the position and attributes for node n0 of edge ei
         */
        vertex<VertexAttribute> v0(size_t ei){ return _nodes[_edges[ei].n0()]; }

	    /**
         * @brief Returns the vertex representing the "second" node in an edge
         * @param ei ID of the edge
         * @return a vertex<VertexAttribute> providing the position and attributes for node n1 of edge ei
         */
        vertex<VertexAttribute> v1(size_t ei){ return _nodes[_edges[ei].n1()]; }

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
        size_t add_edge(size_t inode0, size_t inode1, fiber<VertexAttribute> f, EdgeAttribute a = {}, float epsilon = 0) {

            if (epsilon >= 0) {
                f.remove_duplicates();
            }
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
        void graph_edge(const size_t id, vertex<VertexAttribute>& v0, vertex<VertexAttribute>& v1) {
            const size_t n0 = _edges[id].n0();
            const size_t n1 = _edges[id].n1();

            v0 = _nodes[n0];
            v1 = _nodes[n1];
        }

        /**
         * @brief Returns the path representing the edge as a fiber
         * @param ei index of the edge path to be returned
         * @param include_node_points flag indicating whether or not the nodes will be included in the path (default = true)
         * @return a fiber representing the path of the edge connecting its nodes
         */
        fiber<VertexAttribute> fiber_edge(const size_t ei, bool include_node_points = true) const {
            fiber<VertexAttribute> f = _edges[ei];
            if (include_node_points) {
                f.insert(f.begin(), _nodes[_edges[ei].inodes[0]]);
                f.push_back(_nodes[_edges[ei].inodes[1]]);
            }
            return _edges[ei];
        }

        /**
         * @brief Returns the spatial length of an edge by following the fiber between both nodes
         *
         * @param ei id of the edge to be measured
         * @return the length of the specified edge
         */
        float length(size_t ei) {
            float l = _edges[ei].length(v0(ei), v1(ei));             // calculate the length of the primary fiber
            return l;
        }

        /**
         * @brief Applies a Gaussian smoothing operation to all edges in the network
         * @param sigma standard deviation of the Gaussian kernel
         * @return a new fibernet with all of the fibers smoothed by the kernel
         */
        fibernet smooth_gaussian(float sigma) {
            fibernet smoothed;                                  // create a new fiber network to store the smoothed fibers
            smoothed._nodes = _nodes;                           // store all nodes (their positions will be unchanged)
            for (size_t ei=0; ei<_edges.size(); ei++) {
                edge smoothed_edge = _edges[ei].smooth_gaussian(sigma, _nodes[_edges[ei].n0()], _nodes[_edges[ei].n1()]);
                smoothed._edges.push_back(smoothed_edge);
            }
            return smoothed;
        }
	};
}
