#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <tira/parser.h>

#include <glm/glm.hpp>

namespace tira {

struct triangle {
    glm::vec3 v[3];
    glm::vec3 n;
};

class simplemesh {

    std::vector<triangle> _t;

    glm::vec3 _calculate_normal(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v1v2 = v1 - v2;
        glm::vec3 n = glm::cross(v1v2, v0v1);
        return glm::normalize(n);
    }

    glm::vec3 _load_vertex(tira::parser& p, size_t i) {
        glm::vec3 v;
        v.x = p.get<float>("v", i, 0);
        v.y = p.get<float>("v", i, 1);
        v.z = p.get<float>("v", i, 2);
        return v;
    }

    glm::vec3 _centroid() {
        size_t T = _t.size();           // number of triangles
        glm::vec3 c(0, 0, 0);
        for(size_t ti = 0; ti < T; ti++) {
            for(size_t vi = 0; vi < 3; vi++) {
                c += _t[ti].v[vi];
            }
        }

        float npts = T * 3;
        glm::vec3 result = c / npts;
        return result;
    }

    float _radius(glm::vec3 centroid) {
        size_t T = _t.size();           // number of triangles
        float radius = 0;
        for(size_t ti = 0; ti < T; ti++) {
            for(size_t vi = 0; vi < 3; vi++) {
                glm::vec3 p = _t[ti].v[vi];
                float r = glm::length(p - centroid);
                if(r > radius) radius = r;
            }
        }
        return radius;
    }

    void _scale(float s) {

        const size_t T = _t.size();
        for(size_t ti = 0; ti < T; ti++) {
            for(size_t vi = 0; vi < 3; vi++) {
                _t[ti].v[vi] = s * (_t[ti].v[vi]);
            }
        }
    }

public:

    void load(std::string filename, float scale = 1.0, glm::vec3 position = glm::vec3(0, 0, 0)) {
        tira::parser p(filename);


        std::vector< std::vector<unsigned int> > Faces;
        Faces = p.get<unsigned int>("f");

        std::vector< std::vector<float> > Vertices;
        Vertices = p.get<float>("v");
        size_t F = Faces.size();                // get the number of faces
        _t.resize(F);

        for (size_t fi = 0; fi < F; fi++) {
            size_t v0i = Faces[fi][0] - 1;     // get the vertex IDs for each vertex
            size_t v1i = Faces[fi][1] - 1;          // (remember that OBJ uses 1-based indexing)
            size_t v2i = Faces[fi][2] - 1;

            triangle t;
            t.v[0] = glm::vec3(Vertices[v0i][0], Vertices[v0i][1], Vertices[v0i][2]);
            t.v[1] = glm::vec3(Vertices[v1i][0], Vertices[v1i][1], Vertices[v1i][2]);
            t.v[2] = glm::vec3(Vertices[v2i][0], Vertices[v2i][1], Vertices[v2i][2]);
            t.n = _calculate_normal(t.v[0], t.v[1], t.v[2]);
            _t[fi] = t;
        }
        _scale(scale);
    }

    void boundingsphere(glm::vec3& c, float& r) {
        c = _centroid();
        r = _radius(c);
    }

    // return the number of triangles
    size_t count() {
        return _t.size();
    }

    triangle operator[](size_t ti) {
        return _t[ti];
    }
};





}