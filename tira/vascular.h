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

    struct VesselAttributes {
        float length ;
        float tortuosity;
        float avg_radius;
    };

    typedef fibernet<float, VesselAttributes, unsigned int> vesselnet;
    typedef fiber<float> vessel;

class vascular : public vesselnet {


protected:

    float _length_range[2];

    static constexpr uint8_t major = 1;
    static constexpr uint8_t minor = 0;

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
        return _nodes.size() * sizeof(uint32_t) + _edges.size() * (2 * sizeof(uint32_t) + 3 * sizeof(uint64_t));
    }

    size_t _skeleton_bytes() const {
        size_t node_vertices = _nodes.size() * 4 * 4;       // each node vertex holds four floats (4*4 bytes)
        size_t edge_vertices = _edges.size() * 8;           // each edge stores the number of points in a uint64
        for (size_t ei = 0; ei < _edges.size(); ei++) {
            edge_vertices += _edges[ei].size() * 4 * 4;     // each edge vertex holds four floats (4*4)
        }
        return node_vertices + edge_vertices;               // return the sum of both groups of points
    }
    size_t _surface_bytes() const { return 0; }
    size_t _volume_bytes() const { return 0; }

    void _update_length_range(float length) {
        if (length < _length_range[0])
            _length_range[0] = length;
        if (length > _length_range[1])
            _length_range[1] = length;
    }

    void _calculate_length_range() {
        _length_range[0] = std::numeric_limits<float>::infinity();
        _length_range[1] = 0;
        for (size_t ei=0; ei<_edges.size(); ei++) {
            _update_length_range(_edges[ei].ea().length);
        }
    }

    void _calculate_edge_attributes(size_t ei) {
        _edges[ei].ea().length = length(ei);            // calculate and store the length of the vessel
        _update_length_range(_edges[ei].ea().length);
    }

    void _calculate_attributes() {
        for (size_t ei=0; ei<_edges.size(); ei++) {
            _calculate_edge_attributes(ei);
        }
    }

public:

    void init() {
        _length_range[0] = std::numeric_limits<float>::infinity();
        _length_range[1] = 0;
    }

    vascular(fibernet f) : fibernet(f) {
        init();
        _calculate_attributes();
    }

    vascular() : fibernet() {}

    void save(std::string filename) {
        std::ofstream out(filename, std::ios::binary);                       // create an input file stream
        if (!out.is_open())                                                      // make sure that the file is loaded
            throw std::runtime_error("Could not open file " + filename);    // otherwise throw an exception


        // save the vasc file header
        _Header h;                                                       // create a header structure and fill it with the necessary data
        h.major = major;                                           // get the major and minor version numbers
        h.minor = minor;
        h.num_edges = _edges.size();                                    // number of edges and nodes in the graph
        h.num_nodes = _nodes.size();
        h.skel_offset = sizeof(_Header) + _graph_bytes();                // calculate the size for each of the offsets
        h.surf_offset = h.skel_offset + _skeleton_bytes();
        h.vol_offset = h.surf_offset + _surface_bytes();
        h.write(out);

        // write graph edges
        for (size_t ei = 0; ei < _edges.size(); ei++) {
            size_t n0 = _edges[ei].n0();
            size_t n1 = _edges[ei].n1();
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
        for (size_t ni = 0; ni < _nodes.size(); ni++) {         // iterate through each node
            glm::vec3 p = _nodes[ni];
            float r = _nodes[ni].va();                          // get the vertex attribute (radius)
            out.write((char*)&p[0], 12);
            out.write((char*)&r, 4);
        }

        for (size_t ei = 0; ei < _edges.size(); ei++) {
            size_t n_pts = _edges[ei].size();                   // get the number of points in the edge
            out.write((char*)&n_pts, 8);                        // write it to the file
            for (size_t pi = 0; pi < n_pts; pi++) {
                glm::vec3 p = glm::vec3(_edges[ei][pi]);        // get the point position
                float r = _edges[ei][pi].va();                  // get the radius (vertex attribute)
                out.write((char*)&p[0], 12);                    // write both to the file
                out.write((char*)&r, 4);
            }
        }        
        out.close();
    }

    void load(std::string filename) {
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
            fibernet::edge new_edge(n0, n1);
            _edges.push_back(new_edge);

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
            vertex<float> v(p, r);
            fibernet::node new_node(v, 0);
            _nodes.push_back(new_node);
        }

        for (size_t ei = 0; ei < _edges.size(); ei++) {
            size_t n_pts;                               // get the number of points in the edge
            in.read((char*)&n_pts, 8);                  // write it to the file
            for (size_t pi = 0; pi < n_pts; pi++) {
                glm::vec3 p;                            // get the point position
                float r;                                // get the radius (vertex attribute)
                in.read((char*)&p[0], 12);              // write both to the file
                in.read((char*)&r, 4);
                _edges[ei].push_back(p, r);
            }
        }
        _calculate_attributes();                        // calculate attributes for the vascular network
        in.close();
    }

    /**
     * Add a node to the vascular graph given the nodes medial axis position
     * @param pt is a vec4 structure specifying the (x, y, z) coordinates of the medial axis and the radius at the node
     * @return the internal ID of the added node in the graph
     */
    inline size_t add_node(const glm::vec4 pt) {

        glm::vec3 coord = pt;
        float radius = pt[3];
        vertex<float> new_vertex(coord, radius);

        vesselnet::add_node(new_vertex, 0);
        return vesselnet::_nodes.size() - 1;
    }


    size_t add_edge(size_t inode0, size_t inode1, std::vector<glm::vec4> pts) {

        //create a new fiber to represent the edge
        fiber<float> new_fiber;
        for (size_t pi = 0; pi < pts.size(); pi++) {
            new_fiber.push_back(pts[pi], pts[pi][3]);
        }
        size_t ei = vesselnet::add_edge(inode0, inode1, new_fiber);
        _calculate_edge_attributes(ei);
        return ei;
    }

    /**
     * 
     * @return the number of edges (vessels) in the vasc structure
     */
    size_t edges() const { return _edges.size(); }

    /**
     * 
     * @return the number of nodes (end points and bifurcations) in the vasc structure
     */
    size_t nodes() const { return _nodes.size(); }

    /**
     * @brief Retrieve the vertices representing the nodes of an edge
     *
     * @param id index for the edge that will be retrieved
     * @param v0 first vertex (will be filled with the vertex stored at the first edge node)
     * @param v1 second vertex (will be filled with the vertex stored at the second edge node)
     */
    void edge(size_t id, vertex<float>& v0, vertex<float>& v1) {
        size_t n0 = _edges[id].n0();
        size_t n1 = _edges[id].n1();

        v0 = _nodes[n0];
        v1 = _nodes[n1];
    }

    fiber<> centerline(const size_t id) const {
        return _edges[id];
    }

    vascular smooth(float sigma) {
        fibernet new_fibernet = fibernet::smooth(sigma);
        vascular new_vascular(new_fibernet);
        return new_vascular;
    }
};

}

