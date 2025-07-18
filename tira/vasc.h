﻿#include <vector>
#include <numeric>
#include <string>
#include <fstream>
#include <cstdint>

#include <glm/glm.hpp>


namespace tira {

class vasc {

    static constexpr uint8_t major = 1;
    static constexpr uint8_t minor = 0;

    typedef glm::vec4 medial_pt;

    struct edge {
        std::vector<uint32_t> ipts;      // stores the starting and ending IDs for the points representing the centerline for this edge
        uint32_t inodes[2];              // node indices
    };

    struct node {
        uint32_t ipt;                   // index to the point on the skeleton corresponding to the node
        std::vector<uint32_t> iedges;    // indices of connected edges
    };

    // data structure stores header information from a vasc file
    struct header {
        //uint8_t version[2];            // major and minor version numbers
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

    struct VascEdgeHeader {
        uint32_t node0, node1;
        uint64_t skel_offset;
    };
protected:
    std::vector<medial_pt> _points;
    std::vector<node> _nodes;
    std::vector<edge> _edges;

    size_t _graph_bytes() const;                  // returns the size (in bytes) for each sub-section of the vasc file
    size_t _skeleton_bytes() const;
    size_t _surface_bytes() const;
    size_t _volume_bytes() const;

public:
    void load(const std::string& filename);
    void save(const std::string& filename);
    void obj(const std::string& filename);
    void smooth_paths();
    //void load_vdb(openvdb::FloatGrid::Ptr skeletonGrid);

    // functions used to edit the vasc structure
    uint32_t add_node(glm::vec4 pt);
    uint32_t add_edge(uint32_t inode0, uint32_t inode1, std::vector<glm::vec4> pts);
};

    inline size_t vasc::_graph_bytes() const {
        return _nodes.size() * sizeof(uint32_t) + _edges.size() * (2 * sizeof(uint32_t) + 3 * sizeof(uint64_t));
    }

    inline size_t vasc::_skeleton_bytes() const {
        size_t s = sizeof(uint32_t) + _points.size() * 4 * sizeof(float);
        for (size_t ei = 0; ei < _edges.size(); ei++) {
            s += sizeof(uint32_t) + _edges[ei].ipts.size();
        }
        return s;
    }
    inline size_t vasc::_surface_bytes() const { return 0; }        // Not implemented yet
    inline size_t vasc::_volume_bytes() const { return 0; }         // Not implemented yet

    /**
     * This function loads a .vasc file into the current vasc object
     * @param filename is the name of the .vasc file to load
     */
    void vasc::load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);                       // create an input file stream
        if (!in.is_open())                                                      // make sure that the file is loaded
            throw std::runtime_error("Could not open file " + filename);    // otherwise throw an exception


        //read VASC file header
        header h;
        h.read(in);
        //in.read(reinterpret_cast<char*>(&h), sizeof(header));                    // reads the header directly from the file stream

        _nodes.resize(h.num_nodes);                                             // allocate an array to store node data
        _edges.resize(h.num_edges);                                             // allocate space to store the edge data

        // load the medial axis point ID associated with each node
        for (size_t ni = 0; ni < h.num_nodes; ni++)                                 // iterate through each node
            in.read(reinterpret_cast<char*>(&_nodes[ni].ipt), sizeof(uint32_t));   // load the corresponding ID from the file


        // load the edges (capillaries) from the file
        for (size_t ei = 0; ei < h.num_edges; ei++) {
            in.read(reinterpret_cast<char*>(&_edges[ei].inodes[0]), sizeof(uint32_t));
            in.read(reinterpret_cast<char*>(&_edges[ei].inodes[1]), sizeof(uint32_t));

            // when loading the entire file, these variables are unnecessary (they'll be used
            // later for streaming)
            uint64_t skeleton_offset;
            in.read(reinterpret_cast<char*>(&skeleton_offset), sizeof(uint64_t));
            uint64_t surface_offset;
            in.read(reinterpret_cast<char*>(&surface_offset), sizeof(uint64_t));
            uint64_t volume_offset;
            in.read(reinterpret_cast<char*>(&volume_offset), sizeof(uint64_t));
        }

        // We now have all of the necessary information to generate a graph, so this is a good time
        //  to assign edges to each node structure.
        for (size_t ei = 0; ei < h.num_edges; ei++) {
            _nodes[_edges[ei].inodes[0]].iedges.push_back(ei);
            _nodes[_edges[ei].inodes[1]].iedges.push_back(ei);
        }

        // Load the skeleton data. This consists of a 32-bit integer providing the number of points
        // and then each point stored as a 4d vertex (x, y, z) and radius

        // get the number of medial axis points stored in the file
        uint32_t num_pts;
        in.read(reinterpret_cast<char*>(&num_pts), sizeof(uint32_t));

        _points.resize(num_pts);                                // pre-allocate the array that will store all of the medial axis points
        for (size_t pi = 0; pi < num_pts; pi++) {               // for each point
            in.read(reinterpret_cast<char*>(&_points[pi]), sizeof(medial_pt));      // read the point coordinates and radius into the array
        }

        // Load the indices representing the medial axis for each edge
        for (size_t ei = 0; ei < h.num_edges; ei++) {                           // for each edge
            uint32_t n_pts;
            in.read(reinterpret_cast<char*>(&n_pts), sizeof(uint32_t));       // get the number of medial axis points in the current edge
            _edges[ei].ipts.resize(n_pts);                                      // allocate the length of the medial axis index array for the current edge

            for (size_t epi = 0; epi < n_pts; epi++)                            // for each point in the edge, read the corresponding index
                in.read(reinterpret_cast<char*>(&_edges[ei].ipts[epi]), sizeof(uint32_t));
        }

        in.close();
    }

    /**
     * Saves a vasc file to disk
     * @param filename is the filename
     */
    void vasc::save(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);                       // create an input file stream
        if (!out.is_open())                                                      // make sure that the file is loaded
            throw std::runtime_error("Could not open file " + filename);    // otherwise throw an exception


        // save the vasc file header
        header h;                                                       // create a header structure and fill it with the necessary data
        h.major = major;                                           // get the major and minor version numbers
        h.minor = minor;
        h.num_edges = _edges.size();                                    // number of edges and nodes in the graph
        h.num_nodes = _nodes.size();
        h.skel_offset = sizeof(header) + _graph_bytes();                // calculate the size for each of the offsets
        h.surf_offset = h.skel_offset + _skeleton_bytes();
        h.vol_offset = h.surf_offset + _surface_bytes();
        h.write(out);
        //out.write(reinterpret_cast<char*>(&h), sizeof(header));       // write the header structure to the file


        // write the medial axis point ID associated with each node
        for (size_t ni = 0; ni < _nodes.size(); ni++)                                 // iterate through each node
            out.write(reinterpret_cast<char*>(&_nodes[ni].ipt), sizeof(uint32_t));   // load the corresponding ID from the file


        // write the edges (capillaries) to the file
        size_t skel_edge_offset = h.skel_offset + sizeof(uint32_t) + _points.size() * sizeof(medial_pt);
        size_t surf_edge_offset = h.surf_offset;
        size_t vol_edge_offset = h.vol_offset;

        for (size_t ei = 0; ei < h.num_edges; ei++) {
            out.write(reinterpret_cast<char*>(&_edges[ei].inodes[0]), sizeof(uint32_t));
            out.write(reinterpret_cast<char*>(&_edges[ei].inodes[1]), sizeof(uint32_t));

            // write the data offsets for the edge
            out.write(reinterpret_cast<char*>(&skel_edge_offset), sizeof(uint64_t));
            out.write(reinterpret_cast<char*>(&surf_edge_offset), sizeof(uint64_t));
            out.write(reinterpret_cast<char*>(&vol_edge_offset), sizeof(uint64_t));

            // update the edge offsets for each type of data
            skel_edge_offset += sizeof(uint32_t) + _edges[ei].ipts.size() * sizeof(uint32_t);
            surf_edge_offset += 0;
            vol_edge_offset += 0;
        }

        // Write the skeleton data. This consists of a 32-bit integer providing the number of points
        // and then each point stored as a 4d vertex (x, y, z) and radius

        // get the number of medial axis points stored in the file
        uint32_t num_pts = _points.size();
        out.write(reinterpret_cast<char*>(&num_pts), sizeof(uint32_t));
        out.write(reinterpret_cast<char*>(&_points[0]), sizeof(medial_pt) * num_pts);       // the entire point array can be written at once


        // Write the indices representing the medial axis for each edge
        for (size_t ei = 0; ei < h.num_edges; ei++) {                           // for each edge
            uint32_t n_pts = _edges[ei].ipts.size();
            out.write(reinterpret_cast<char*>(&n_pts), sizeof(uint32_t));       // get the number of medial axis points in the current edge

            out.write(reinterpret_cast<char*>(&_edges[ei].ipts[0]), sizeof(uint32_t));  // write the index array to disk
        }

        // In the future we'll have to put the surface and volume writing code here

        out.close();
    }

    /**
     * Add a node to the vascular graph given the nodes medial axis position
     * @param pt is a vec4 structure specifying the (x, y, z) coordinates of the medial axis and the radius at the node
     * @return the internal ID of the added node in the graph
     */
    uint32_t vasc::add_node(glm::vec4 pt) {
        _points.push_back(pt);                          // add the new node point to the list of medial axis points
        node new_node;                                  // create a new node structure
        new_node.ipt = _points.size() - 1;              // add the ID for the new medial point to the node structure
        _nodes.push_back(new_node);                     // push the node structure to the nodes array

        return _nodes.size() - 1;
    }

    /**
     * Add an edge to the graph given its connected nodes and a set of points on the medial axis
     * @param inode0 the first node in the graph (associated with the medial axis point nearest the first point in pts)
     * @param inode1 the second node in the graph (associated with the medial axis point closest to the last point in pts)
     * @param pts vector of vec4 medial axis points consisting of an (x, y, z) coordinate and radius
     * @return the internal ID of the edge in the graph
     */
    uint32_t vasc::add_edge(uint32_t inode0, uint32_t inode1, std::vector<glm::vec4> pts) {

        edge new_edge;                                  // create a new edge structure
        new_edge.inodes[0] = inode0;                    // set the two node IDs based on user input
        new_edge.inodes[1] = inode1;
        new_edge.ipts.resize(pts.size());               // resize the medial axis point index array to match the number of points in the edge
        std::iota(new_edge.ipts.begin(), new_edge.ipts.end(), _points.size());  // create a sequential list of IDs starting with the last index in _points
        _edges.push_back(new_edge);

        _points.insert(_points.end(), pts.begin(), pts.end());      // insert all of the new points for this edge into the _points array

        return _edges.size() - 1;
    }


}

