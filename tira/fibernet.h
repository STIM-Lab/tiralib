#pragma once

#include "fiber.h"
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

namespace tira {

    /**
     * @brief      This edge class extends a fiber to include a user-defined edge attribute
     * and a list of connected edges
     */
    template<typename VertexAttributeType, typename EdgeAttributeType>
    class edge : public fiber<VertexAttributeType> {
    protected:
        EdgeAttributeType m_edge_attribute;
        size_t m_node_indices[2];              // node indices

    public:

        /**
         * @brief      Constructs a new edge given an existing fiber, two nodes, and an edge attribute
         *
         * @param[in]  f       existing fiber containing a series of geometric vertices and their attributes
         * @param[in]  n0      first node connected by this edge (closest to the first vertex in the fiber)
         * @param[in]  n1      second node connected to this edge (closest to the last vertex in the fiber)
         * @param[in]  attrib  user-defined attribute to store information within each edge
         */
        edge(fiber<VertexAttributeType> f, size_t n0, size_t n1, EdgeAttributeType attrib = {}) : fiber<VertexAttributeType>(f) {
            m_node_indices[0] = n0;
            m_node_indices[1] = n1;
            m_edge_attribute = attrib;
        }

        /**
         * @brief      Creates an empty edge linking two nodes
         *
         * @param[in]  n0      first node ID connected by this edge
         * @param[in]  n1      second node ID connected by this edge
         * @param[in]  attrib  user-defined attribute that stores information at this edge
         */
        edge(size_t n0, size_t n1, EdgeAttributeType attrib = {}) {
            m_node_indices[0] = n0;
            m_node_indices[1] = n1;
            m_edge_attribute = attrib;
        }

        /**
         * @brief       Creates a new edge using an existing edge and new fiber - Used to update the geometry of the edge.
         *
         * @param f     New fiber geometry
         * @param e     Existing edge, where the graph and attribute parameters will be copied
         */
        edge(fiber<VertexAttributeType> f, edge& e) : fiber<VertexAttributeType>(f) {
            m_node_indices[0] = e.m_node_indices[0];
            m_node_indices[1] = e.m_node_indices[1];
            m_edge_attribute = e._ea;
        }

        /**
         * @brief      Assigns or adjusts nodes connected by this edge.
         *
         * @param[in]  n0    sets the first node connected by this edge
         * @param[in]  n1    sets the second node connected by this edge
         */
        void nodes(size_t n0, size_t n1) {
            m_node_indices[0] = n0;
            m_node_indices[1] = n1;
        }

        /**
         * @brief      Retrieve the first node ID
         *
         * @return     index into the _nodes vector
         */
        size_t NodeIndex0() const { return m_node_indices[0]; }

        /**
         * @brief      Retrieve the first node ID
         *
         * @return     index into the _nodes vector
         */
        size_t NodeIndex1() const { return m_node_indices[1]; }

        /**
         * @brief      Returns the user-specified edge attribute
         *
         * @return     The edge attribute.
         */
        EdgeAttributeType& EdgeAttribute() { return m_edge_attribute; }

        /**
         * @brief      Assigns an attribute to this edge
         *
         * @param[in]  attrib  The attribute
         */
        void EdgeAttribute(EdgeAttributeType attrib) { m_edge_attribute = attrib; }

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
        edge Smooth(float sigma, vertex<VertexAttributeType> node_v0, vertex<VertexAttributeType> node_v1) {
            fiber<VertexAttributeType> original_fiber = *this;

            original_fiber.insert(original_fiber.begin(), node_v0);
            original_fiber.push_back(node_v1);
            fiber<VertexAttributeType> smoothed_fiber = original_fiber.Smooth(sigma);
            smoothed_fiber.erase(smoothed_fiber.begin());
            smoothed_fiber.pop_back();

            edge smoothed_edge(smoothed_fiber, m_node_indices[0], m_node_indices[1]);
            return smoothed_edge;
        }

        float Length(vertex<VertexAttributeType> v0, vertex<VertexAttributeType> v1) const {
            fiber<VertexAttributeType> original_fiber = *this;
            original_fiber.insert(original_fiber.begin(), v0);
            original_fiber.push_back(v1);
            return original_fiber.Length();
        }

    };

    /**
     * @brief      Defines a class for a node, which connects multiple fibers and consists of a geometric point and attributes
     */
    template<typename VertexAttributeType, typename NodeAttributeType>
    class node : public vertex<VertexAttributeType> {
    protected:
        NodeAttributeType m_node_attribute;
        std::vector<size_t> m_edge_indices;   // indices of connected edges

    public:

        node() : vertex<VertexAttributeType>() {}

        node(const node& n) : vertex<VertexAttributeType>(n) {
            m_node_attribute = n.m_node_attribute;
            m_edge_indices = n.m_edge_indices;
        }

        node& operator=(node x) {
            vertex<VertexAttributeType>::operator=(x);      // call the vertex assignment operator
            m_node_attribute = x.m_node_attribute;
            m_edge_indices = x.m_edge_indices;
            return *this;
        }

        /**
         * @brief      Create a new node from an existing vertex and a node-specific attribute
         *
         * @param[in]  v       Existing vertex
         * @param[in]  attrib  Node attribute
         */
        node(vertex<VertexAttributeType> v, NodeAttributeType attrib) : vertex<VertexAttributeType>(v) {
            m_node_attribute = attrib;
        }



        /**
         * @brief      Connects this node to an edge. This should be done in parity with connecting edges to nodes
         *
         * @param[in]  ei    ID of the edge connecting to this node
         */
        void AddEdgeIndex(size_t ei) { m_edge_indices.push_back(ei); }

        /**
         * @brief      Retrieve the node-specific attribute
         *
         * @return     The node attribute.
         */
        NodeAttributeType NodeAttribute() { return m_node_attribute; }

        /**
         * @brief      Assign or update an attribute for this node
         *
         * @param[in]  na    The new value for the node attribute
         */
        void NodeAttribute(NodeAttributeType na) { m_node_attribute = na; }

        void ClearEdgeIndices() { m_edge_indices.clear(); }

        ///////////////////////////////////////////////

        /**
         * @brief returns the number of edges connected to this node (which is  the degree).
        */
        size_t Degree() const {
            return m_edge_indices.size(); // where degegree = number of connected edges
        }


        /**
         * @brief removes an edge index from the node�s list of connected edges.
         *        this is called when an edge is deleted from the graph.
        */

        void RemoveEdge(size_t ei) {
            // remove all occurrences of edge index ei from the edge list
            m_edge_indices.erase(
                std::remove(m_edge_indices.begin(), m_edge_indices.end(), ei),
                m_edge_indices.end()
            );
        }

        /**
         * Returns the list of edges attached to this node
         * @return indices for the attached edges
         */
        std::vector<size_t> Edges() const { return m_edge_indices; }

        /////////////////////////////////

    };

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
	template <typename VertexAttributeType = float, 
        typename EdgeAttributeType = unsigned int,
        typename NodeAttributeType = unsigned int>
	class fibernet {
	public:
	    typedef node<VertexAttributeType, NodeAttributeType> fnode;
	    typedef edge<VertexAttributeType, EdgeAttributeType> fedge;

    protected:
        /**
         * List of nodes in this network, essentially modeled as a graph
         */
        std::vector<fnode> m_nodes;

        /**
         * List of edges in this network, modeled as a mathematical graph
         */
        std::vector<fedge> m_edges;


        /**
         * Assign connected edges to each node
         */
        void m_UpdateNodeEdgeList() {

            for (size_t ni = 0; ni < m_nodes.size(); ni++)                 // clear the edge indices from all nodes
                m_nodes[ni].ClearEdgeIndices();

            for (size_t ei = 0; ei < m_edges.size(); ei++) {        // add edge indices to each associated node
                size_t ni0 = m_edges[ei].NodeIndex0();
                size_t ni1 = m_edges[ei].NodeIndex1();
                m_nodes[ni0].AddEdgeIndex(ei);
                m_nodes[ni1].AddEdgeIndex(ei);
            }
	    }

    public:
        /**
         * @brief Returns the vertex representing the "first" node in an edge
         * @param ei ID of the edge
         * @return a vertex<VertexAttribute> providing the position and attributes for node n0 of edge ei
         */
        vertex<VertexAttributeType> NodeVertex0(size_t ei){ return m_nodes[m_edges[ei].NodeIndex0()]; }

	    /**
         * @brief Returns the vertex representing the "second" node in an edge
         * @param ei ID of the edge
         * @return a vertex<VertexAttribute> providing the position and attributes for node n1 of edge ei
         */
        vertex<VertexAttributeType> NodeVertex1(size_t ei){ return m_nodes[m_edges[ei].NodeIndex1()]; }

        /**
         * @brief      Adds a node to the network using a geometric position and attribute
         *
         * @param[in]  v     vertex describing the geometric position of the node
         * @param[in]  n     node-specific attributes
         *
         * @return     { description_of_the_return_value }
         */
        size_t AddNode(vertex<VertexAttributeType> v, NodeAttributeType n) {
            node new_node(v, n);
            m_nodes.push_back(new_node);
            return m_nodes.size() - 1;
        }

	    void BoundingBox(size_t edge_id, glm::vec3& aabb_min, glm::vec3& aabb_max) {

            size_t n0i = m_edges[edge_id].NodeIndex0();
            size_t n1i = m_edges[edge_id].NodeIndex1();

            bool result = m_edges[edge_id].BoundingBox(aabb_min, aabb_max);

            if (!result) aabb_min = m_nodes[n0i];

            aabb_min = glm::min(m_nodes[n0i], aabb_min);
            aabb_max = glm::max(m_nodes[n0i], aabb_max);

            aabb_min = glm::min(m_nodes[n1i], aabb_min);
            aabb_max = glm::max(m_nodes[n1i], aabb_max);

        }

        /**
         * Add an edge to the graph given its connected nodes and fiber geometry
         * @param inode0 the first node in the graph (associated with the medial axis point nearest the first point in pts)
         * @param inode1 the second node in the graph (associated with the medial axis point closest to the last point in pts)
         * @param f is the fiber geometry
         * @return the internal ID of the edge in the graph
         */
        size_t AddEdge(size_t inode0, size_t inode1, fiber<VertexAttributeType> f, EdgeAttributeType a = {}) {

            f.RemoveDuplicates();
            
            edge new_edge(f, inode0, inode1, a);                                  // create a new edge structure
            m_edges.push_back(new_edge);
            
            size_t idx = m_edges.size() - 1;
            m_nodes[inode0].AddEdgeIndex(idx);
            m_nodes[inode1].AddEdgeIndex(idx);

            return m_edges.size() - 1;
        }

         /**
         * Returns the 3D coordinates and radii associated with both nodes associated with an edge.
         * @param id identifier for the edge that will be returned
         * @return std::pair storing the coordinates and radii of both nodes connected by the edge
         */
        void Endpoints(const size_t id, vertex<VertexAttributeType>& v0, vertex<VertexAttributeType>& v1) {
            const size_t n0 = m_edges[id].NodeIndex0();
            const size_t n1 = m_edges[id].NodeIndex1();

            v0 = m_nodes[n0];
            v1 = m_nodes[n1];
        }

        /**
         * @brief Returns the path representing the edge as a fiber
         * @param ei index of the edge path to be returned
         * @param include_node_points flag indicating whether or not the nodes will be included in the path (default = true)
         * @return a fiber representing the path of the edge connecting its nodes
         */
        fiber<VertexAttributeType> Fiber(const size_t ei, bool include_node_points = true) const {
            fiber<VertexAttributeType> f = m_edges[ei];
            if (include_node_points) {
                f.insert(f.begin(), m_nodes[m_edges[ei].inodes[0]]);
                f.push_back(m_nodes[m_edges[ei].inodes[1]]);
            }
            return m_edges[ei];
        }

        /**
         * @brief Returns the spatial length of an edge by following the fiber between both nodes
         *
         * @param ei id of the edge to be measured
         * @return the length of the specified edge
         */
        float Length(size_t ei) {
            float l = m_edges[ei].Length(NodeVertex0(ei), NodeVertex1(ei));             // calculate the length of the primary fiber
            return l;
        }

	    float ChordLength(size_t ei) {
            vertex<VertexAttributeType> v0 = m_nodes[m_edges[ei].NodeIndex0()];
            vertex<VertexAttributeType> v1 = m_nodes[m_edges[ei].NodeIndex1()];

            return glm::length(v1 - v0);
        }

        /**
         * @brief Applies a Gaussian smoothing operation to all edges in the network
         * @param sigma standard deviation of the Gaussian kernel
         * @return a new fibernet with all of the fibers smoothed by the kernel
         */
        fibernet Smooth(float sigma) {
            fibernet smoothed;                                  // create a new fiber network to store the smoothed fibers
            smoothed.m_nodes = m_nodes;                           // store all nodes (their positions will be unchanged)
            for (size_t ei=0; ei< m_edges.size(); ei++) {
                edge smoothed_edge = m_edges[ei].Smooth(sigma, m_nodes[m_edges[ei].NodeIndex0()], m_nodes[m_edges[ei].NodeIndex1()]);
                smoothed.m_edges.push_back(smoothed_edge);
            }
            return smoothed;
        }

        ////////////////////Delete and merge implementation ////////////////////////////

        /**
         * @brief deletes a node if its degree is zero
         * shifts all node indices in edges above this node
         *
         * @param node_id: index of the node to delete
         * return: true if deleted, false otherwise
        */

        bool DeleteNode(size_t node_id) {

            if (node_id >= m_nodes.size()) return false;

            if (m_nodes[node_id].Degree() > 0) throw std::runtime_error("cannot delete node with nonzero degree");

            for (auto& e : m_edges) { // shift edge node references that come after the deleted node

                if (e.NodeIndex0() > node_id) e.nodes(e.NodeIndex0() - 1, e.NodeIndex1());

                if (e.NodeIndex1() > node_id) e.nodes(e.NodeIndex0(), e.NodeIndex1() - 1);
            }
            m_nodes.erase(m_nodes.begin() + node_id); // erase the node
            return true;
        }

        const std::vector<fnode>& InternalNodes() const { return m_nodes; }
        const std::vector<fedge>& InternalEdges() const { return m_edges; }

        /**
         * @brief deletes an edge from the graph and updates all references.
         * if the connected nodes become isolated, they are deleted
         *
         * @param edge_id: index of the edge to delete
        */

        void DeleteEdge(size_t edge_id) {

            if (edge_id >= m_edges.size()) return;

            size_t n0 = m_edges[edge_id].NodeIndex0(), n1 = m_edges[edge_id].NodeIndex1();

            m_nodes[n0].RemoveEdge(edge_id); m_nodes[n1].RemoveEdge(edge_id); // remove edge reference from both nodes

            for (auto& node : m_nodes) {
                std::vector<size_t> updated;

                for (auto ei : node.Edges()) {

                    if (ei > edge_id) updated.push_back(ei - 1);

                    else if (ei < edge_id) updated.push_back(ei);
                    // else skip (ei == edge_id, already removed)
                }
                node.ClearEdgeIndices();

                for (auto ei : updated)
                    node.AddEdgeIndex(ei);
            }


            m_edges.erase(m_edges.begin() + edge_id); // remove the edge

            if (m_nodes[n0].Degree() == 0) DeleteNode(n0); // delete n0 if isolated

            if (m_nodes[n1].Degree() == 0 && n1 != n0) DeleteNode(n1); // delete n1 if isolated and not same as n0
        }


        /**
        * @brief merges two edges that share a single node (which must have degree 2)
        * deletes both edges, deletes shared node (if disconnected), and adds new merged edge
        *
        * @param e1: index of the first edge
        * param e2: index of the second edge
       */

        /*
        
         Start at a degree-2 node that hasn't been visited yet.

         expand forward through the graph following only other degree-2 nodes.

         then expand backward the same way.

         now have a list of connected nodes  and their connecting edges.

         if this path has 3 or more nodes, can merge it:

         concatenate their fiber geometry in order.

         Delete all  degree 2 ones between the ends.

         delete all edges.
  
         add a new edge from the first to the last node with the merged fiber.

         Repeat for all degree 2 chains.
        
        */
        void MergeDegree2Chains() {

            std::vector<bool> visited_node(m_nodes.size(), false);
            std::vector<bool> visited_edge(m_edges.size(), false);

            auto is_mergeable = [&](size_t node) {
                return node < m_nodes.size() && m_nodes[node].Degree() == 2;
                };

            size_t totalMerged = 0;

            for (size_t ni = 0; ni < m_nodes.size(); ++ni) {
                if (!is_mergeable(ni) || visited_node[ni])
                    continue;

                std::vector<size_t> chain_nodes{ ni };
                std::vector<size_t> chain_edges;

                // Expand forward
                size_t curr = ni;
                size_t prev = static_cast<size_t>(-1);
                while (true) {
                    visited_node[curr] = true;
                    const auto& edges = m_nodes[curr].Edges();
                    size_t next = static_cast<size_t>(-1), next_edge = static_cast<size_t>(-1);

                    for (size_t ei : edges) {
                        if (visited_edge[ei]) continue;

                        const auto& edge = m_edges[ei];
                        size_t n0 = edge.NodeIndex0(), n1 = edge.NodeIndex1();
                        size_t other = (n0 == curr) ? n1 : n0;

                        if (other != prev && is_mergeable(other)) {
                            next = other;
                            next_edge = ei;
                            break;
                        }
                    }

     
                    if (next_edge == static_cast<size_t>(-1)) break;

                    chain_edges.push_back(next_edge);
                    curr = next;
                    prev = chain_nodes.back();
                    chain_nodes.push_back(curr);
                }

                // Expand backward
                curr = ni;
                prev = static_cast<size_t>(-1);
                while (true) {
                    const auto& edges = m_nodes[curr].Edges();
                    size_t next = static_cast<size_t>(-1), next_edge = static_cast<size_t>(-1);

                    for (size_t ei : edges) {
                        if (visited_edge[ei]) continue;

                        const auto& edge = m_edges[ei];
                        size_t n0 = edge.NodeIndex0(), n1 = edge.NodeIndex1();
                        size_t other = (n0 == curr) ? n1 : n0;

                        if (other != prev && is_mergeable(other)) {
                            next = other;
                            next_edge = ei;
                            break;
                        }
                    }

                    if (next_edge == static_cast<size_t>(-1)) break;

                    chain_edges.insert(chain_edges.begin(), next_edge);
                    curr = next;
                    prev = chain_nodes.front();
                    chain_nodes.insert(chain_nodes.begin(), curr);
                    visited_node[curr] = true;
                }

                // Validate
                if (chain_nodes.size() < 3) continue;

                size_t start_node = chain_nodes.front();
                size_t end_node = chain_nodes.back();

                // Build merged fiber
                tira::fiber<float> merged_fiber;
                for (size_t i = 0; i < chain_edges.size(); ++i) {
                    const auto& edge = m_edges[chain_edges[i]];
                    fiber<float> segment = edge;

                    if (edge.NodeIndex1() == chain_nodes[i]) {
                        std::reverse(segment.begin(), segment.end());
                    }

                    if (merged_fiber.empty())
                        merged_fiber.insert(merged_fiber.end(), segment.begin(), segment.end());
                    else
                        merged_fiber.insert(merged_fiber.end(), segment.begin() + 1, segment.end()); // avoid duplication
                }

                // Delete edges (reverse order)
                std::sort(chain_edges.begin(), chain_edges.end(), std::greater<size_t>());
                for (size_t ei : chain_edges) {
                    visited_edge[ei] = true;
                    DeleteEdge(ei);
                }

                // Delete internal nodes
                for (size_t i = 1; i + 1 < chain_nodes.size(); ++i) {
                    size_t nid = chain_nodes[i];
                    if (m_nodes[nid].Degree() == 0) {
                        DeleteNode(nid);
                    }
                }

                AddEdge(start_node, end_node, merged_fiber);
                totalMerged++;
            }

            std::cout << "Merged " << totalMerged << " degree 2 chains." << std::endl;
        }

        /////////////////////////////////////////////////////////////////

        /**
         * @brief Removes an edge from the edge list and updates connected node indices accordingly.
         *
         * @param edge_index Index of the edge to remove.
        */

        void RemoveEdge(size_t edge_index) {
            if (edge_index >= m_edges.size()) return;

            const auto& edge = m_edges[edge_index];

            // decrease the degree of the connected nodes
            m_nodes[edge.NodeIndex0()].RemoveEdge(edge_index);
            m_nodes[edge.NodeIndex1()].RemoveEdge(edge_index);

            // erase the edge from the list.......
            m_edges.erase(m_edges.begin() + edge_index);

        }

        /**
         * @brief Prints the number of nodes with the specified degree.
         *
         * @param degree the target node degree to search for.
        */

        void QueryDegree(int degree) const {
            size_t count = 0;

            for (const auto& node : m_nodes) {
                if (node.Degree() == degree) {
                    count++;
                }
            }

            if (count == 0) {
                std::cout << "no nodes with degree " << degree << "." << std::endl;
            }
            else {
                std::cout << count << " nodes with degree " << degree << "." << std::endl;
            }
        }

        /**
         * @brief check if the edge at index `edge_idx` is a spine
         * @param edge_idx Index of the edge
         * @param min_len minimum length threshold
         * @param max_len maximum length threshold
         * @return True if the edge is a spine, false otherwise
        */

        bool IsSpine(size_t edge_idx) const {

            const auto& edge = m_edges[edge_idx];

            /*float len = edge.Length(m_nodes[edge.NodeIndex0()], m_nodes[edge.NodeIndex1()]);

            if (len < min_len || len > max_len)
                return false;                                                             // reject if length is outside allowed range
            */
            // get the two endpoint nodes
            const auto& n0 = m_nodes[edge.NodeIndex0()];
            const auto& n1 = m_nodes[edge.NodeIndex1()];

            return (n0.Degree() == 1 || n1.Degree() == 1);                                // check if either node has degree 1
        }


        /**
         * @brief Selects edges that are "spines" based on length and degree.
         *        Supports AND/OR logic with previous selection and optional removal.
         *
         * @param current   Optional vector of already selected edge indices
         * @param and_op        If true, perform AND (intersection); else OR (union)
         * @return          vector of edge indices connected to at least one degree-1 node
        */

        std::vector<size_t> QuerySpines(const std::vector<size_t>& current = {}, bool and_op = false)
        {
            std::vector<size_t> result;

            if (and_op) {
                for (size_t ei : current) {
                    if (IsSpine(ei)) {
                        result.push_back(ei);
                    }
                }
            }
            else {
                for (size_t ei = 0; ei < m_edges.size(); ++ei) {
                    if (IsSpine(ei)) {
                        result.push_back(ei);
                    }
                }

                // OR: add original current
                result.insert(result.end(), current.begin(), current.end());

                // Remove duplicates
                std::sort(result.begin(), result.end());
                result.erase(std::unique(result.begin(), result.end()), result.end());
            }

            /*if (remove) {
                std::sort(result.begin(), result.end(), std::greater<size_t>());
                for (size_t ei : result) {
                    this->RemoveEdge(ei);
                }
            }*/

            return result;
        }


        ///////////////////////////////////////////////////

        /**
         * @brief returns a new fibernet with the edges in edge_indices removed,
         * preserving the original network.
         * @param edge_indices list of edge indices to delete.
         * @return a new fibernet with the specified edges removed.
         */
        fibernet DeleteEdges(const std::vector<size_t>& edge_indices) const {
            std::vector<bool> to_delete(m_edges.size(), false);
            for (size_t ei : edge_indices) {
                if (ei < m_edges.size()) to_delete[ei] = true;
            }
            fibernet new_net;
            new_net.m_nodes = m_nodes; // copy all nodes

            // copy only edges that are not marked for deletion
            for (size_t ei = 0; ei < m_edges.size(); ++ei) {
                if (!to_delete[ei]) {
                    new_net.m_edges.push_back(m_edges[ei]);
                }
            }

            // rebuild node edge indices
            for (auto& node : new_net.m_nodes)
                node.ClearEdgeIndices();

            for (size_t ei = 0; ei < new_net.m_edges.size(); ++ei) {
                new_net.m_nodes[new_net.m_edges[ei].NodeIndex0()].AddEdgeIndex(ei);
                new_net.m_nodes[new_net.m_edges[ei].NodeIndex1()].AddEdgeIndex(ei);
            }
            return new_net;
        }


        /**
         * @brief Selects all edges whose total fiber length falls within the given range [low, high].
         *        Can be chained with previous results using AND (intersection) or OR (union).
         * @param low The minimum allowed length for an edge .
         * @param high The maximum allowed length for an edge .
         * @param current Optional vector of previously selected edge indices.
         * @param op      If true: AND (restrict to edges in 'current' AND in range);
         *                If false: OR (include any edge in 'current' OR in range).
         * @return        A vector of edge indices that satisfy the length criteria.
        */
        std::vector<size_t> QueryLength( float low, float high, const std::vector<size_t>& current = {}, bool op = false) const {
            std::vector<size_t> result;
            std::vector<bool> already_in(m_edges.size(), false);

            if (op && !current.empty()) {
                // AND: only look at the current selection
                for (size_t i : current) {
                    if (i < m_edges.size()) {
                        float len = m_edges[i].Length(
                            m_nodes[m_edges[i].NodeIndex0()],
                            m_nodes[m_edges[i].NodeIndex1()]);
                        if (len >= low && len <= high)
                            result.push_back(i);
                    }
                }
            }
            else {
                // OR: go over all edges, add anything that matches or is already in current
                for (size_t i : current)
                    if (i < m_edges.size())
                        already_in[i] = true;

                for (size_t ei = 0; ei < m_edges.size(); ++ei) {
                    float len = m_edges[ei].Length(
                        m_nodes[m_edges[ei].NodeIndex0()],
                        m_nodes[m_edges[ei].NodeIndex1()]);
                    if ((len >= low && len <= high) && !already_in[ei]) {
                        result.push_back(ei);
                    }
                    else if (already_in[ei]) {
                        result.push_back(ei);
                    }

                }
            }
            return result;
        }


	    /**
         * @brief Calculate the mean absolute curvature of an edge
         * @param edge_idx index of the edge to be analyzed
         * @return mean absolute curvature of the edge
        */
	    float MeanCurvature(size_t edge_idx) {

            edge<VertexAttributeType, EdgeAttributeType>& current_edge = m_edges[edge_idx];                                // get the fiber (fiber<float>)
            if (current_edge.size() < 3) return 0.0f;                                 // need at least 3 points for second derivative

            std::vector<float> kappa = current_edge.Curvature();                     // call the curvature function

            float sum_abs_curvature = 0.0f;
            for (float k : kappa)
                sum_abs_curvature += std::abs(k);                            // accumulate absolute curvature

            return sum_abs_curvature / static_cast<float>(kappa.size());     // compute mean
        }


        /**
         * @brief selects all edges whose mean absolute curvature (tortuosity) falls within the range [kmin, kmax].
         *         with previous results using AND (intersection) or OR (union).
         * @param kmin    The minimum allowed tortuosity value.
         * @param kmax    The maximum allowed tortuosity value.
         * @param current Optional vector of previously selected edge indices.
         * @param op      If true: AND (restrict to edges in 'current' AND in range);
         *                If false: OR (include any edge in 'current' OR in range).
         * @return        A vector of edge indices that satisfy the tortuosity criteria.
        */
        std::vector<size_t> QueryMeanCurvature(float kmin, float kmax, const std::vector<size_t>& current = {}, bool op = false) {

            std::vector<size_t> result;

            // AND operation
            if (op) {
                for (size_t ci = 0; ci < current.size(); ci++) {
                    float tort = MeanCurvature(current[ci]);
                    if (tort >= kmin && tort <= kmax) {
                        result.push_back(current[ci]);
                    }
                }
                return result;
            }

            // OR operation
            for (size_t ei = 0; ei < m_edges.size(); ei++) {
                float tort = MeanCurvature(ei);
                if (tort >= kmin && tort <= kmax) {
                    result.push_back(ei);
                }
            }

            // combine the result vector with the input vector
            result.insert(result.end(), current.begin(), current.end());

            //remove duplicates
            std::sort(result.begin(), result.end());
            std::unique(result.begin(), result.end());


            return result;
        }

        /**
         * @brief Arc/Chord ratio for an edge . Uses polyline length (arc) over straight-line (chord).
        */
        float ArcChordRatio(size_t edge_idx) const {

            const auto& e = m_edges[edge_idx];
            const auto& v0 = m_nodes[e.NodeIndex0()];
            const auto& v1 = m_nodes[e.NodeIndex1()];

            float arc = e.Length(v0, v1);                                // polyline length including interior points
            float chord = glm::length(glm::vec3(v1) - glm::vec3(v0));      // straight line

            if (chord <= 1e-6f) return std::numeric_limits<float>::infinity(); // avoid divide-by-zero
            return arc / chord;
        }

        /**
         * @brief selects all edges whose tortuosity falls within the given range [tmin, tmax].
         *        with previous results using AND (intersection) or OR (union).
         * @param tmin    The minimum allowed tortuosity for an edge.
         * @param tmax    The maximum allowed tortuosity for an edge.
         * @param current Optional vector of previously selected edge indices.
         * @param op      If true: AND (restrict to edges in 'current' AND in range);
         *                If false: OR (include any edge in 'current' OR in range).
         * @return        A vector of edge indices that satisfy the tortuosity criteria.
        */

        std::vector<size_t> QueryVesselArcChord( float tmin, float tmax, const std::vector<size_t>& current = {}, bool op = false)
        {
            std::vector<size_t> result;

            // AND: restrict to 'current'
            if (op) {
                for (size_t ei : current) {
                    float r = ArcChordRatio(ei);
                    if (r >= tmin && r <= tmax)
                        result.push_back(ei);
                }
                return result;
            }

            // OR: scan all edges
            for (size_t ei = 0; ei < m_edges.size(); ++ei) {
                float r = ArcChordRatio(ei);
                if (r >= tmin && r <= tmax)
                    result.push_back(ei);
            }

            // combine the result vector with the input vector
            result.insert(result.end(), current.begin(), current.end());

            //remove duplicates
            std::sort(result.begin(), result.end());
            result.erase(std::unique(result.begin(), result.end()), result.end());

            return result;
        }

	    size_t Degree(size_t node_id) {
            return m_nodes[node_id].Degree();
        }

	    size_t Node0(size_t edge_id) {
            return m_edges[edge_id].NodeIndex0();
        }
	    size_t Node1(size_t edge_id) {
            return m_edges[edge_id].NodeIndex1();
        }

	    /**
	     * Returns a histogram of the degrees of all nodes
         * @return
         */
        std::vector<size_t> CountDegrees() {

            std::vector<size_t> histogram;

            for (size_t ni = 0; ni < m_nodes.size(); ni++) {
                size_t degree = m_nodes[ni].Degree();
                if (degree >= histogram.size())
                    histogram.resize(degree + 1);
                histogram[degree]++;
            }
            return histogram;
        }

        /**
         * Retrieve all node indices connected to the specified edges
         * @param edges list of edge indices connected to the desired nodes
         * @return unique node indices connected to the specified edges (duplicates are removed)
         */
        std::vector<size_t> Nodes(const std::vector<size_t>& edges) {

            std::vector<size_t> connected_nodes;
            connected_nodes.reserve(edges.size() * 2);          // reserve two nodes for each edge

            for (size_t i = 0; i < edges.size(); i++) {
                size_t ei = edges[i];                           // get an edge index

                connected_nodes.emplace_back(m_edges[ei].NodeIndex0());     // put both nodes connected to this edge in the list
                connected_nodes.emplace_back(m_edges[ei].NodeIndex1());
            }

            // remove duplicates
            std::sort(connected_nodes.begin(), connected_nodes.end());
            connected_nodes.erase( std::unique(connected_nodes.begin(), connected_nodes.end()), connected_nodes.end() );
            return connected_nodes;
        }

	    std::vector<size_t> Edges(const std::vector<size_t>& nodes) {
            std::vector<size_t> connected_edges;
            connected_edges.reserve(nodes.size() * 3);      // reserve three edges for each node (seems reasonable for vasculature)

            for (size_t i = 0; i < nodes.size(); i++) {
                size_t ni = nodes[i];

                std::vector<size_t> attached_edges = m_nodes[ni].Edges();
                connected_edges.insert(connected_edges.end(), attached_edges.begin(), attached_edges.end());
            }

            // remove duplicates
            std::sort(connected_edges.begin(), connected_edges.end());
            connected_edges.erase( std::unique(connected_edges.begin(), connected_edges.end()), connected_edges.end() );
            return connected_edges;
        }

	};
}
