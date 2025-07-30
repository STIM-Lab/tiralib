#pragma once
#include <vector>
#include <numeric>
#include <string>
#include <fstream>
#include <cstdint>
#include <limits>

#include <glm/glm.hpp>

#include "fibernet.h"


namespace tira {

    struct VesselAttributesType {
        float length ;
        float tortuosity;
        float avg_radius;
    };

    typedef fibernet<float, VesselAttributesType, unsigned int> vesselnet;
    typedef fiber<float> vessel;

class vascular : public vesselnet {


protected:

    float m_length_range[2];

    static constexpr uint8_t s_major = 1;
    static constexpr uint8_t s_minor = 0;

    // data structure stores header information from a vasc file
    struct _Header {
        uint8_t major;
        uint8_t minor;
        uint32_t num_nodes;             // number of nodes in the vascular network
        uint32_t num_edges;             // total number of edges in the vascular network (corresponds to microvessesl)
        uint64_t skel_offset;           // number of bytes from the beginning of the file where the skeleton information starts
        uint64_t surf_offset;           // number of bytes from the beginning of the file where the surface information starts
        uint64_t vol_offset;           // number of bytes from the beginning of the file where the volume information starts

        void write(std::ofstream& fout) {
            fout.write(reinterpret_cast<char*>(&major), sizeof(uint8_t));
            fout.write(reinterpret_cast<char*>(&minor), sizeof(uint8_t));
            fout.write(reinterpret_cast<char*>(&num_nodes), sizeof(uint32_t));
            fout.write(reinterpret_cast<char*>(&num_edges), sizeof(uint32_t));
            fout.write(reinterpret_cast<char*>(&skel_offset), sizeof(uint64_t));
            fout.write(reinterpret_cast<char*>(&surf_offset), sizeof(uint64_t));
            fout.write(reinterpret_cast<char*>(&vol_offset), sizeof(uint64_t));
        }

        void read(std::ifstream& fin) {
            fin.read(reinterpret_cast<char*>(&major), sizeof(uint8_t));
            fin.read(reinterpret_cast<char*>(&minor), sizeof(uint8_t));
            fin.read(reinterpret_cast<char*>(&num_nodes), sizeof(uint32_t));
            fin.read(reinterpret_cast<char*>(&num_edges), sizeof(uint32_t));
            fin.read(reinterpret_cast<char*>(&skel_offset), sizeof(uint64_t));
            fin.read(reinterpret_cast<char*>(&surf_offset), sizeof(uint64_t));
            fin.read(reinterpret_cast<char*>(&vol_offset), sizeof(uint64_t));
        }
    };

    struct _VascEdgeHeader {
        uint32_t node0, node1;
        uint64_t skel_offset;
    };

    size_t _graph_bytes() const {                  // returns the size (in bytes) for each sub-section of the vasc file
        return m_nodes.size() * sizeof(uint32_t) + m_edges.size() * (2 * sizeof(uint32_t) + 3 * sizeof(uint64_t));
    }

    size_t _skeleton_bytes() const {
        size_t node_vertices = m_nodes.size() * 4 * 4;       // each node vertex holds four floats (4*4 bytes)
        size_t edge_vertices = m_edges.size() * 8;           // each edge stores the number of points in a uint64
        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            edge_vertices += m_edges[ei].size() * 4 * 4;     // each edge vertex holds four floats (4*4)
        }
        return node_vertices + edge_vertices;               // return the sum of both groups of points
    }
    size_t _surface_bytes() const { return 0; }
    size_t _volume_bytes() const { return 0; }

    void _update_length_range(float length) {
        if (length < m_length_range[0])
            m_length_range[0] = length;
        if (length > m_length_range[1])
            m_length_range[1] = length;
    }

    void _calculate_length_range() {
        m_length_range[0] = std::numeric_limits<float>::infinity();
        m_length_range[1] = 0;
        for (size_t ei=0; ei< m_edges.size(); ei++) {
            _update_length_range(m_edges[ei].ea().length);
        }
    }

    void _calculate_edge_attributes(size_t ei) {
        m_edges[ei].ea().length = Length(ei);            // calculate and store the length of the vessel
        _update_length_range(m_edges[ei].ea().length);
    }

    void _calculate_attributes() {
        for (size_t ei=0; ei< m_edges.size(); ei++) {
            _calculate_edge_attributes(ei);
        }
    }

    /**
     * Removes duplicated vertices within fibers and between fibers and nodes. These shouldn't occur, so a warning
     * is output whenever it happens.
     */
    void _remove_duplicate_points() {
        for (size_t ei=0; ei< m_edges.size(); ei++) {
            size_t removed = m_edges[ei].RemoveDuplicates();
            if (removed != 0) {
                std::cout<<"tira::vascular WARNING: "<<removed<<" duplicate points detected in the fiber associated with edge "<<ei<<std::endl;
            }

            // test for duplicates between the edge fiber and its nodes
            vertex<float> v0 = m_nodes[m_edges[ei].NodeIndex0()];
            if (v0 == m_edges[ei][0]) {
                m_edges[ei].erase(m_edges[ei].begin());
                std::cout<<"tira::vascular WARNING: node 0 is duplicated with the fiber in edge "<<ei<<std::endl;
            }
            vertex<float> v1 = m_nodes[m_edges[ei].NodeIndex1()];
            if (v1 == m_edges[ei].back()) {
                m_edges[ei].pop_back();
                std::cout<<"tira::vascular WARNING: node 1 is duplicated with the fiber in edge "<<ei<<std::endl;
            }
        }
    }

public:

    void Init() {
        m_length_range[0] = std::numeric_limits<float>::infinity();
        m_length_range[1] = 0;
    }

    vascular(fibernet f) : fibernet(f) {
        Init();
        _calculate_attributes();
    }

    vascular() : fibernet() {}

    void Save(std::string filename) {
        std::ofstream out(filename, std::ios::binary);                       // create an input file stream
        if (!out.is_open())                                                      // make sure that the file is loaded
            throw std::runtime_error("Could not open file " + filename);    // otherwise throw an exception


        // save the vasc file header
        _Header h;                                                       // create a header structure and fill it with the necessary data
        h.major = s_major;                                           // get the major and minor version numbers
        h.minor = s_minor;
        h.num_edges = m_edges.size();                                    // number of edges and nodes in the graph
        h.num_nodes = m_nodes.size();
        h.skel_offset = sizeof(_Header) + _graph_bytes();                // calculate the size for each of the offsets
        h.surf_offset = h.skel_offset + _skeleton_bytes();
        h.vol_offset = h.surf_offset + _surface_bytes();
        h.write(out);

        // write graph edges
        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            size_t n0 = m_edges[ei].NodeIndex0();
            size_t n1 = m_edges[ei].NodeIndex1();
            out.write((char*)&n0, 8);
            out.write((char*)&n1, 8);
            size_t skel_offset = 0;
            size_t surf_offset = 0;
            size_t vol_offset = 0;
            out.write((char*)&skel_offset, 8);
            out.write((char*)&surf_offset, 8);
            out.write((char*)&vol_offset, 8);
        }

        // write skeleton
        for (size_t ni = 0; ni < m_nodes.size(); ni++) {         // iterate through each node
            glm::vec3 p = m_nodes[ni];
            float r = m_nodes[ni].Attribute();                          // get the vertex attribute (radius)
            out.write((char*)&p[0], 12);
            out.write((char*)&r, 4);
        }

        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            size_t n_pts = m_edges[ei].size();                   // get the number of points in the edge
            out.write((char*)&n_pts, 8);                        // write it to the file
            for (size_t pi = 0; pi < n_pts; pi++) {
                glm::vec3 p = glm::vec3(m_edges[ei][pi]);        // get the point position
                float r = m_edges[ei][pi].Attribute();                  // get the radius (vertex attribute)
                out.write((char*)&p[0], 12);                    // write both to the file
                out.write((char*)&r, 4);
            }
        }        
        out.close();
    }

    void Load(std::string filename) {
        std::ifstream in(filename, std::ios::binary);                       // create an input file stream
        if (!in.is_open())                                                      // make sure that the file is loaded
            throw std::runtime_error("Could not open file " + filename);    // otherwise throw an exception


        // load the vasc file header
        _Header h;                                                       // create a header structure and fill it with the necessary data
        h.read(in);

        // read graph edges
        for (size_t ei = 0; ei < h.num_edges; ei++) {
                        
            size_t n0, n1;
            in.read((char*)&n0, 8);
            in.read((char*)&n1, 8);
            fibernet::fedge new_edge(n0, n1);
            m_edges.push_back(new_edge);

            size_t skel_offset;
            size_t surf_offset;
            size_t vol_offset;
            in.read((char*)&skel_offset, 8);
            in.read((char*)&surf_offset, 8);
            in.read((char*)&vol_offset, 8);
        }

        // read skeleton
        for (size_t ni = 0; ni < h.num_nodes; ni++) {           // iterate through each node
            glm::vec3 p;
            float r;                                            // get the vertex attribute (radius)
            in.read((char*)&p[0], 12);
            in.read((char*)&r, 4);
            vertex v(p, r);
            fnode new_node(v, 0);
            m_nodes.push_back(new_node);
        }

        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            size_t n_pts;                               // get the number of points in the edge
            in.read((char*)&n_pts, 8);                  // write it to the file
            for (size_t pi = 0; pi < n_pts; pi++) {
                glm::vec3 p;                            // get the point position
                float r;                                // get the radius (vertex attribute)
                in.read((char*)&p[0], 12);              // read both from the file
                in.read((char*)&r, 4);
                m_edges[ei].AddLastVertex(p, r);
            }
        }
        _remove_duplicate_points();
        _calculate_attributes();                        // calculate attributes for the vascular network
        in.close();
    }

    /**
     * Add a node to the vascular graph given the nodes medial axis position
     * @param pt is a vec4 structure specifying the (x, y, z) coordinates of the medial axis and the radius at the node
     * @return the internal ID of the added node in the graph
     */
    inline size_t AddNode(const glm::vec4 pt) {

        glm::vec3 coord = pt;
        float radius = pt[3];
        vertex<float> new_vertex(coord, radius);

        vesselnet::AddNode(new_vertex, 0);
        return vesselnet::m_nodes.size() - 1;
    }


    size_t AddEdge(size_t inode0, size_t inode1, std::vector<glm::vec4> pts) {

        //create a new fiber to represent the edge
        fiber<float> new_fiber;
        for (size_t pi = 0; pi < pts.size(); pi++) {
            new_fiber.AddLastVertex(pts[pi], pts[pi][3]);
        }
        size_t ei = vesselnet::AddEdge(inode0, inode1, new_fiber);
        _calculate_edge_attributes(ei);
        return ei;
    }

    /**
     * 
     * @return the number of edges (vessels) in the vasc structure
     */
    size_t NumEdges() const { return m_edges.size(); }

    /**
     * 
     * @return the number of nodes (end points and bifurcations) in the vasc structure
     */
    size_t NumNodes() const { return m_nodes.size(); }

    /**
     * @brief returns a new fibernet with the edges in edge_indices removed,
     * preserving the original network.
     * @param edge_indices list of edge indices to delete.
     * @return a new fibernet with the specified edges removed.
     */
    vascular DeleteEdges(const std::vector<size_t>& edge_indices) const {
        std::vector<bool> to_delete(m_edges.size(), false);
        for (size_t ei : edge_indices) {
            if (ei < m_edges.size()) to_delete[ei] = true;
        }
        vascular new_net;
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
    std::vector<size_t> SelectEdgeLength( float low, float high, const std::vector<size_t>& current = {}, bool op = false
    ) const {
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
     * @brief selects all edges whose mean radius falls within the given range [rmin, rmax].
     *        can be chained with previous results using AND (intersection) or OR (union).
     * @param rmin    The minimum allowed mean radius for an edge.
     * @param rmax    The maximum allowed mean radius for an edge.
     * @param current Optional vector of previously selected edge indices.
     * @param op      If true: AND (restrict to edges in 'current' AND in range);
     *                If false: OR (include any edge in 'current' OR in range).
     * @return        A vector of edge indices that satisfy the radius criteria.
    */
    std::vector<size_t> SelectEdgeRadius( float rmin, float rmax, const std::vector<size_t>& current = {}, bool op = false
    ) const {
        std::vector<size_t> result;
        std::vector<bool> already_in(m_edges.size(), false);

        if (op && !current.empty()) {
            // AND: only look at the current selection
            for (size_t i : current) {
                if (i < m_edges.size()) {
                    const auto& edge = m_edges[i];
                    float sum = 0;
                    for (const auto& pt : edge) sum += pt.Attribute();
                    float mean = edge.empty() ? 0.0f : sum / edge.size();
                    if (mean >= rmin && mean <= rmax)
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
                const auto& edge = m_edges[ei];
                float sum = 0.0f;
                for (const auto& pt : edge)
                    sum += pt.Attribute();

                float mean = 0.0f;
                if (!edge.empty()) {
                    mean = sum / static_cast<float>(edge.size());
                }

                if ((mean >= rmin && mean <= rmax) && !already_in[ei]) {
                    result.push_back(ei);
                }
                else if (already_in[ei]) {
                    result.push_back(ei);
                }


            }
        }
        return result;
    }


    fiber<> Centerline(size_t id, bool include_nodes = true) {
        fiber<float> c = m_edges[id];
        if (include_nodes) {
            vertex<float> v0 = m_nodes[vesselnet::m_edges[id].NodeIndex0()];
            c.insert(c.begin(), v0);
            vertex<float> v1 = m_nodes[vesselnet::m_edges[id].NodeIndex1()];
            c.push_back(v1);
        }
        return c;
    }

    vascular Smooth(float sigma) {
        fibernet new_fibernet = fibernet::Smooth(sigma);
        vascular new_vascular(new_fibernet);
        return new_vascular;
    }
};

}

