#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <unordered_set>

struct medial_pt {
    float p[3]; // position: x, y, z
    float r;    // radius
};

struct edge {
    std::vector<uint32_t> pts; // indices into points
    uint32_t nodes[2];          // node indices
};

struct node {
    unsigned pt;                // index into points
    std::vector<uint32_t> edges; // indices of connected edges
};

struct VascHeader {
    uint32_t num_vertices;
    uint32_t num_edges;
    uint64_t skel_offset;
};

struct VascEdgeHeader {
    uint32_t node0, node1;
    uint64_t skel_offset;
};

struct CoordHash {
    std::size_t operator()(const openvdb::Coord& c) const {
        return std::hash<int>()(c.x()) ^ (std::hash<int>()(c.y()) << 1) ^ (std::hash<int>()(c.z()) << 2);
    }
};


class vasc {
protected:
    std::vector<medial_pt> _points;
    std::vector<node> _nodes;
    std::vector<edge> _edges;

public:
    void load(const std::string& filename);
    void save(const std::string& filename);
    void obj(const std::string& filename);
    void smooth_paths();
    void load_vdb(openvdb::FloatGrid::Ptr skeletonGrid);
};

/*
This function is loading an in-memory graph structure
_points = all 3D coordinates with radius
_nodes = key junction points 
_edges = curved lines between nodes, stored as paths made of points

*/
void vasc::load(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "failed to open VASC file: " << filename << "\n";
        return;
    }

    //read VASC file header
    VascHeader header;
    in.read(reinterpret_cast<char*>(&header), sizeof(header)); // contains num_vertices, num_edges, etc.

    // node to point mapping 
    std::vector<uint32_t> nodePts(header.num_vertices);
    for (auto& id : nodePts)
        in.read(reinterpret_cast<char*>(&id), sizeof(uint32_t));  // one uint32 per node // reinterpret_cast<char*> is for reading binary data directly into memory.

    // read edge headers (node pair + skeleton offset) 
    std::vector<VascEdgeHeader> edgeHeaders(header.num_edges);
    in.read(reinterpret_cast<char*>(edgeHeaders.data()), sizeof(VascEdgeHeader) * header.num_edges);

    //  _nodes with point references 
    _nodes.resize(header.num_vertices);
    for (size_t i = 0; i < _nodes.size(); ++i)
        _nodes[i].pt = nodePts[i];

    //reading per_edge skeleton data 
    _edges.resize(header.num_edges);
    for (size_t i = 0; i < edgeHeaders.size(); ++i) {
        // seek to this edge's skeleton data offset
        in.seekg(edgeHeaders[i].skel_offset);

        // how many points are on this edge
        uint32_t num_pts;
        in.read(reinterpret_cast<char*>(&num_pts), sizeof(num_pts));

        // fill edge data
        edge& e = _edges[i];
        e.nodes[0] = edgeHeaders[i].node0;
        e.nodes[1] = edgeHeaders[i].node1;

        for (uint32_t j = 0; j < num_pts; ++j) {
            medial_pt pt;
            in.read(reinterpret_cast<char*>(&pt), sizeof(pt));  // read x, y, z, r
            e.pts.push_back(_points.size());                    // store global point index
            _points.push_back(pt);                              // store point
        }

        // link edge to its two nodes
        _nodes[e.nodes[0]].edges.push_back(i);
        _nodes[e.nodes[1]].edges.push_back(i);
    }

    in.close();
}

/*
It writes in-memory (made of _points, _nodes, _edges) to a .vasc file 
*/

void vasc::save(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "failed to write VASC file: " << filename << "\n";
        return;
    }

    //  initial VASC header (with placeholder skel_offset = 0) 
    VascHeader header{ (uint32_t)_nodes.size(), (uint32_t)_edges.size(), 0 };
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // write node table (each node's point index into _points) 
    for (const node& n : _nodes) {
        out.write(reinterpret_cast<const char*>(&n.pt), sizeof(uint32_t));
    }

    // write placeholder edge headers (node0, node1, skel_offset = 0) 
    std::vector<std::streampos> edge_offset_positions;
    for (const edge& e : _edges) {
        VascEdgeHeader h{ e.nodes[0], e.nodes[1], 0 };
        edge_offset_positions.push_back(out.tellp());  // save where to patch later
        out.write(reinterpret_cast<const char*>(&h), sizeof(h));
    }

    //write skeleton path data (per edge point sequences) 
    header.skel_offset = static_cast<uint64_t>(out.tellp()); // store true offset
    std::vector<uint64_t> skel_offsets;

    for (const edge& e : _edges) {
        skel_offsets.push_back(static_cast<uint64_t>(out.tellp()));

        // write number of points in this edge
        uint32_t npts = (uint32_t)e.pts.size();
        out.write(reinterpret_cast<const char*>(&npts), sizeof(npts));

        // write each medial_pt (x, y, z, radius) from _points
        for (uint32_t pid : e.pts) {
            out.write(reinterpret_cast<const char*>(&_points[pid]), sizeof(medial_pt));
        }
    }

    // backpatch header with correct skel_offset 
    out.seekp(0);
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // backpatch edge headers with correct skel_offsets 
    for (size_t i = 0; i < _edges.size(); ++i) {
        out.seekp(static_cast<std::streamoff>(edge_offset_positions[i]) + 8);
        out.write(reinterpret_cast<const char*>(&skel_offsets[i]), sizeof(uint64_t));
    }

    out.close();
}


void vasc::obj(const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to write OBJ: " << filename << "\n";
        return;
    }
    for (const medial_pt& pt : _points)
        out << "v " << pt.p[0] << " " << pt.p[1] << " " << pt.p[2] << "\n";
    for (const edge& e : _edges) {
        for (size_t i = 1; i < e.pts.size(); ++i)
            out << "l " << e.pts[i - 1] + 1 << " " << e.pts[i] + 1 << "\n";
    }
    out.close();
}

void vasc::smooth_paths() {
    for (edge& e : _edges) {
        if (e.pts.size() < 3) continue;  // Not enough points to smooth

        // create a copy of original points
        std::vector<medial_pt> new_pts(e.pts.size());

        // first and last points stay the same
        new_pts[0] = _points[e.pts[0]];
        new_pts.back() = _points[e.pts.back()];

        // apply 3 point averaging
        for (size_t i = 1; i + 1 < e.pts.size(); ++i) {
            const medial_pt& prev = _points[e.pts[i - 1]];
            const medial_pt& curr = _points[e.pts[i]];
            const medial_pt& next = _points[e.pts[i + 1]];

            medial_pt smoothed;
            smoothed.p[0] = (prev.p[0] + curr.p[0] + next.p[0]) / 3.0f;
            smoothed.p[1] = (prev.p[1] + curr.p[1] + next.p[1]) / 3.0f;
            smoothed.p[2] = (prev.p[2] + curr.p[2] + next.p[2]) / 3.0f;
            smoothed.r = curr.r;  // keep original radius

            new_pts[i] = smoothed;
        }

        // update global _points and e.pts with new values
        for (size_t i = 0; i < e.pts.size(); ++i) {
            _points[e.pts[i]] = new_pts[i];
        }
    }
}

///////////// steps /////////////////////
/*
Extract active voxels with value > 0 from the input VDB grid.

Convert them to world coordinates, assign an index and default radius.

Build an adjacency list by checking 26-connected neighbors.

Detect node points, voxels with degree not equal to 2.

Trace edges between nodes by walking through degree-2 voxels.

Store the graph as _nodes, _edges, and _points.


*/

void vasc::load_vdb(openvdb::FloatGrid::Ptr grid) {
    using Coord = openvdb::Coord;
    using Vec3f = openvdb::Vec3f;

    // map from voxel coordinate to unique index (used to track point IDs)
    std::unordered_map<openvdb::Coord, uint32_t, CoordHash> coordToIndex;

    // store the world space coordinates of active skeleton voxels
    std::vector<Vec3f> worldPoints;

    // extract active voxels from grid 
    for (auto iter = grid->cbeginValueOn(); iter; ++iter) {
        if (iter.getValue() <= 0.0f) continue; // skip background
        Coord ijk = iter.getCoord();
        Vec3f xyz = grid->transform().indexToWorld(ijk); // convert to world coords
        coordToIndex[ijk] = static_cast<uint32_t>(worldPoints.size());
        worldPoints.push_back(xyz);
    }

    // store medial points with default radius 
    for (const auto& p : worldPoints) {
        medial_pt mp;
        mp.p[0] = p.x(); mp.p[1] = p.y(); mp.p[2] = p.z();
        mp.r = 1.0f; // set default radius (can be updated later via SDF)
        _points.push_back(mp);
    }

    // create adjacency graph for all points 
    std::vector<std::vector<int>> adjacency(_points.size());
    std::vector<Coord> neighbors;

    // build 26-neighborhood (6-face + 12-edge + 8-corner)
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
                if (dx || dy || dz) // exclude (0,0,0)
                    neighbors.emplace_back(dx, dy, dz);

    // for each point, record neighbors that are also part of skeleton
    for (const auto& [coord, idx] : coordToIndex) {
        for (const auto& n : neighbors) {
            Coord neighborCoord = coord + n;
            auto it = coordToIndex.find(neighborCoord);
            if (it != coordToIndex.end()) {
                adjacency[idx].push_back(it->second);
            }
        }
    }

    // identify graph nodes (i.e. voxels with degree ≠ 2) 
    std::unordered_map<int, uint32_t> vertexToNodeID;
    std::vector<int> nodeIndices;
    uint32_t nextNodeID = 0;

    for (size_t i = 0; i < adjacency.size(); ++i) {
        int deg = adjacency[i].size();
        if (deg != 2) {
            // this voxel is a node (endpoint, branch point, like this..)
            vertexToNodeID[i] = nextNodeID++;
            node n;
            n.pt = i;
            _nodes.push_back(n);
            nodeIndices.push_back(i);
        }
    }

    // trace edges between nodes 
    std::unordered_set<std::pair<int, int>, std::function<size_t(const std::pair<int, int>&)>> visited(
        0, [](const std::pair<int, int>& p) {
            return std::hash<int>()(p.first) ^ std::hash<int>()(p.second);
        });

    for (int node : nodeIndices) {
        for (int nbr : adjacency[node]) {
            auto key = std::minmax(node, nbr);
            if (visited.count(key)) continue; // already processed this edge

            edge e;
            e.nodes[0] = vertexToNodeID[node]; // start node ID

            std::vector<uint32_t> path;
            path.push_back(node);

            int prev = node;
            int curr = nbr;

            // walk along the line until another node is found
            while (vertexToNodeID.find(curr) == vertexToNodeID.end()) {
                path.push_back(curr);
                const auto& nbrs = adjacency[curr];
                int next = (nbrs[0] == prev) ? nbrs[1] : nbrs[0];
                prev = curr;
                curr = next;
            }

            // reached the destination node
            path.push_back(curr);
            e.nodes[1] = vertexToNodeID[curr]; // end node ID
            e.pts = path;

            // store edge
            uint32_t edgeID = static_cast<uint32_t>(_edges.size());
            _edges.push_back(e);

            // attach edge to both endpoint nodes
            _nodes[e.nodes[0]].edges.push_back(edgeID);
            _nodes[e.nodes[1]].edges.push_back(edgeID);

            visited.insert(key); // mark as visited
        }
    }

    std::cout << "Loaded VDB skeleton into vasc structure.\n";
}



