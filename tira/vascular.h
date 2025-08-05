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
        float length;
        float avg_curvature;
        float avg_radius;
        float arc_cord_ratio;
        float volume;
        float surface_area;
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
            _update_length_range(m_edges[ei].EdgeAttribute().length);
        }
    }

    void _calculate_edge_attributes(size_t ei) {
        m_edges[ei].EdgeAttribute().length = Length(ei);            // calculate and store the length of the vessel
        _update_length_range(m_edges[ei].EdgeAttribute().length);
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

    /**
     * Calculates the length of every vessel in the vascular network and stores that length as a vessel attribute
     * This will allow us to access the length easily any time we need it.
     */
    void m_CalculateLengths() {
        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            vertex<float> v0 = m_nodes[m_edges[ei].NodeIndex0()];
            vertex<float> v1 = m_nodes[m_edges[ei].NodeIndex1()];
            m_edges[ei].EdgeAttribute().length = m_edges[ei].Length(v0, v1);
        }
    }


    /**
     * @brief Calculate the surface area of a vessel using a segment-wise integration approach.
     *        The surface area is estimated by summing 2π·r·Δs for each segment, where the radius is averaged
     *        between adjacent points (including node endpoints).
     * @param edge_idx index of the edge to be analyzed
     * @return estimated lateral surface area of the edge
    */

    float m_CalculateSurfaceArea(size_t edge_idx) {

        float surface_area = 0.0f;                                                    // initialize running surface area
        const auto& edge = m_edges[edge_idx];                                         // reference the edge fiber
        const size_t num_pts = edge.size();                                           // get number of intermediate fiber points

        // initialize the first point (node0)
        glm::vec3 p0 = m_nodes[edge.NodeIndex0()];                                    // get the position of the first node
        float r0 = m_nodes[edge.NodeIndex0()].Attribute();                            // get the radius of the first node

        // loop over all fiber points and node1
        for (size_t pi = 0; pi <= num_pts; ++pi) {

            glm::vec3 p1;                                                              // coordinates of the next point
            float r1;                                                                  // radius of the next point

            if (pi < num_pts) {
                p1 = edge[pi];                                                         // fiber point
                r1 = edge[pi].Attribute();                                             // radius at the fiber point
            }
            else {
                p1 = m_nodes[edge.NodeIndex1()];                                       // final node
                r1 = m_nodes[edge.NodeIndex1()].Attribute();                           // final node radius
            }

            float seg_len = glm::length(p1 - p0);                                      // segment length
            float r_avg = 0.5f * (r0 + r1);                                            // average radius of the segment

            surface_area += 2.0f * static_cast<float>(M_PI) * r_avg * seg_len;        // accumulate lateral surface area

            p0 = p1;                                                                   // update current point
            r0 = r1;                                                                   // update current radius
        }

        return surface_area;                                                           // return the total surface area
    }

    /**
     * Calculate the surface area of every vessel in the network and stores them in the corresponding
     * vessel attribute.
     */
    void m_CalculateSurfaceAreas() {
        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            // calculate the surface area of the edge and store it as an edge attribute
            m_edges[ei].EdgeAttribute().surface_area = m_CalculateSurfaceArea(ei);
        }
    }

    void m_CalculateArcChords() {
        for (size_t ei = 0; ei < m_edges.size(); ei++) {
            // calculate the surface area of the edge and store it as an edge attribute
            float chord_length = ChordLength(ei);
            float arc_length = Length(ei);

            m_edges[ei].EdgeAttribute().arc_cord_ratio = arc_length / chord_length;
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
     * @brief Calculate the average radius of an edge (including node points)
     * @param edge_idx index of the edge to be analyzed
     * @return average radius of the edge
     */
    float AverageRadius(size_t edge_idx) {

        float sum_radii = m_nodes[m_edges[edge_idx].NodeIndex0()].Attribute();      // initialize a running sum with the radius of the first node vertex

        size_t num_pts = m_edges[edge_idx].size();                                  // get the number of vertices in the fiber
        for (size_t pi = 0; pi < num_pts; pi++)                                     // for each vertex
            sum_radii += m_edges[edge_idx][pi].Attribute();                         // add the corresponding radius into the running sum

        sum_radii += m_nodes[m_edges[edge_idx].NodeIndex1()].Attribute();           // add the radius of the second node vertex

        float radius = sum_radii / (num_pts + 2);                                   // calculate the average from the running sum
        return radius;                                                              // return the radius
    }

    /**
     * @brief Calculate the vessel volume of an edge using a segment-wise integration approach.
     *        The volume is estimated by summing the cylindrical volume of each segment ,
     *        where the radius is averaged between adjacent points (including node endpoints).
     * @param edge_idx index of the edge to be analyzed
     * @return estimated vessel volume of the edge
    */
    float VesselVolume(size_t edge_idx) {

        float volume = 0.0f;                                                         // initialize running volume
        const auto& edge = m_edges[edge_idx];                                        // reference the edge fiber
        const size_t num_pts = edge.size();                                          // get number of intermediate points

        // initialize the first point (node0)
        glm::vec3 p0 = m_nodes[edge.NodeIndex0()];                                   // get the position of the first node
        float r0 = m_nodes[edge.NodeIndex0()].Attribute();                           // get the radius of the first node

        // loop over all fiber points and node1
        for (size_t pi = 0; pi <= num_pts; ++pi) {

            glm::vec3 p1;                                                             // coordinates of the next point
            float r1;                                                                 // radius of the next point

            if (pi < num_pts) {
                p1 = edge[pi];                                                        // fiber point
                r1 = edge[pi].Attribute();                                            // radius at the fiber point
            }
            else {
                p1 = m_nodes[edge.NodeIndex1()];                                      // final node
                r1 = m_nodes[edge.NodeIndex1()].Attribute();                          // final node radius
            }

            float seg_len = glm::length(p1 - p0);                                     // segment length
            float r_avg = 0.5f * (r0 + r1);                                           // average radius of the segment

            volume += static_cast<float>(M_PI) * r_avg * r_avg * seg_len;            // add the cylindrical segment volume

            p0 = p1;                                                                  // update current point
            r0 = r1;                                                                  // update current radius
        }

        return volume;                                                                // return the final accumulated volume
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
    std::vector<size_t> QueryVesselRadius(float rmin, float rmax, const std::vector<size_t>& current = {}, bool op = false) {

        std::vector<size_t> result;                                 // initialize a vector to store the result

        // AND operation
        if (op) {                                                   // an AND operation just requires looking at the edges in the "current" vector
            for (size_t ci = 0; ci < current.size(); ci++) {        // for each edge in the current vector
                float c_radius = AverageRadius(current[ci]);  // calculate the radius
                if (c_radius >= rmin && c_radius <= rmax) {         // check to see if it's within the specified range
                    result.push_back(current[ci]);                  // if so, push it into the result vector
                }
            }
            return result;                                          // return the result vector - we're done
        }

        // OR operation
        for (size_t ei = 0; ei < m_edges.size(); ei++) {            // an OR operation requires looking at every edge in the network
            float e_radius = AverageRadius(ei);               // calculate the average radius of the edge
            if (e_radius >= rmin && e_radius <= rmax) {             // if it's within the specified range
                result.push_back(ei);                               // add it to the result vector
            }
        }

        // combine the result vector with the input vector
        result.insert(result.end(), current.begin(), current.end());

        // remove duplicates
        std::sort(result.begin(), result.end());
        std::unique(result.begin(), result.end());


        return result;
    }

    /**
     * @brief selects all edges whose vessel volume falls within the given range [vmin, vmax].
     *         with previous results using AND (intersection) or OR (union).
     * @param vmin    The minimum allowed vessel volume for an edge.
     * @param vmax    The maximum allowed vessel volume for an edge.
     * @param current Optional vector of previously selected edge indices.
     * @param op      If true: AND (restrict to edges in 'current' AND in range);
     *                If false: OR (include any edge in 'current' OR in range).
     * @return        A vector of edge indices that satisfy the volume criteria.
    */

    std::vector<size_t> QueryVesselVolume(float vmin, float vmax, const std::vector<size_t>& current = {}, bool op = false) {

        std::vector<size_t> result;                                     // initialize a vector to store the result

        // AND operation
        if (op) {                                                       // an AND operation just requires looking at the edges in the "current" vector
            for (size_t ci = 0; ci < current.size(); ci++) {            // for each edge in the current vector
                float c_volume = VesselVolume(current[ci]);             // calculate the volume
                if (c_volume >= vmin && c_volume <= vmax) {             // check to see if it's within the specified range
                    result.push_back(current[ci]);                      // if so, push it into the result vector
                }
            }
            return result;                                              // return the result vector - we're done
        }

        // OR operation
        for (size_t ei = 0; ei < m_edges.size(); ei++) {                // an OR operation requires looking at every edge in the network
            float e_volume = VesselVolume(ei);                          // calculate the volume of the edge
            if (e_volume >= vmin && e_volume <= vmax) {                 // if it's within the specified range
                result.push_back(ei);                                   // add it to the result vector
            }
        }

        // combine the result vector with the input vector
        result.insert(result.end(), current.begin(), current.end());

        // remove duplicates
        std::sort(result.begin(), result.end());
        std::unique(result.begin(), result.end());

        return result;                                                  // return the list of matching edge indices
    }




    /**
     * @brief selects all edges whose surface area falls within the given range [amin, amax].
     *        can be chained with previous results using AND (intersection) or OR (union).
     * @param amin    The minimum allowed surface area for an edge.
     * @param amax    The maximum allowed surface area for an edge.
     * @param current Optional vector of previously selected edge indices.
     * @param op      If true: AND (restrict to edges in 'current' AND in range);
     *                If false: OR (include any edge in 'current' OR in range).
     * @return        A vector of edge indices that satisfy the surface area criteria.
     */

    std::vector<size_t> QueryVesselSurfaceArea(float amin, float amax, const std::vector<size_t>& current = {}, bool op = false) {

        std::vector<size_t> result;                                      // initialize a vector to store the result

        // AND operation
        if (op) {                                                        // an AND operation just requires looking at the edges in the "current" vector
            for (size_t ci = 0; ci < current.size(); ci++) {             // for each edge in the current vector
                //float c_area = VesselSurfaceArea(current[ci]);           // calculate the surface area
                float c_area = m_edges[current[ci]].EdgeAttribute().surface_area;
                if (c_area >= amin && c_area <= amax) {                  // check to see if it's within the specified range
                    result.push_back(current[ci]);                       // if so, push it into the result vector
                }
            }
            return result;                                               // return the result vector - we're done
        }

        // OR operation
        for (size_t ei = 0; ei < m_edges.size(); ei++) {                 // an OR operation requires looking at every edge in the network
            float e_area = m_edges[current[ei]].EdgeAttribute().surface_area;                        // calculate the surface area of the edge
            if (e_area >= amin && e_area <= amax) {                      // if it's within the specified range
                result.push_back(ei);                                    // add it to the result vector
            }
        }

        // combine the result vector with the input vector
        result.insert(result.end(), current.begin(), current.end());

        // remove duplicates
        std::sort(result.begin(), result.end());
        std::unique(result.begin(), result.end());

        return result;                                                   // return the list of matching edge indices
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
    std::vector<size_t> QueryVesselArcChord(float tmin, float tmax, const std::vector<size_t>& current = {}, bool op = false) {

        std::vector<size_t> result;                                     // initialize a vector to store the result

        // AND operation
        if (op) {                                                       // an AND operation just requires looking at the edges in the "current" vector
            for (size_t ci = 0; ci < current.size(); ci++) {            // for each edge in the current vector
                float c_tort = m_edges[current[ci]].EdgeAttribute().arc_cord_ratio;                 // calculate the tortuosity
                if (c_tort >= tmin && c_tort <= tmax) {                 // check to see if it's within the specified range
                    result.push_back(current[ci]);                      // if so, push it into the result vector
                }
            }
            return result;                                              // return the result vector - we're done
        }

        // OR operation
        for (size_t ei = 0; ei < m_edges.size(); ei++) {                // an OR operation requires looking at every edge in the network
            float e_tort = m_edges[current[ei]].EdgeAttribute().arc_cord_ratio;                              // calculate the tortuosity of the edge
            if (e_tort >= tmin && e_tort <= tmax) {                     // if it's within the specified range
                result.push_back(ei);                                   // add it to the result vector
            }
        }

        // combine the result vector with the input vector
        result.insert(result.end(), current.begin(), current.end());

        // remove duplicates
        std::sort(result.begin(), result.end());
        std::unique(result.begin(), result.end());

        return result;                                                  // return the list of matching edge indices
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

    /**
     * Calculate all of the network attributes and store them for fast access
     */
    void CalculateAttributes() {
        m_CalculateSurfaceAreas();                  // calculate the surface area for each vessel
        m_CalculateLengths();                       // calculate the length of each vessel
    }

    vascular Smooth(float sigma) {
        fibernet new_fibernet = fibernet::Smooth(sigma);
        vascular new_vascular(new_fibernet);
        return new_vascular;
    }
};

}

