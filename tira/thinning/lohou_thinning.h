#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>

#include <tira/volume.h>

// input arguments
std::string in_inputname;
std::string in_outputname;
std::string in_templates;
int in_cuda;
unsigned int in_foreground;

struct coord {
    int dx;
    int dy;
    int dz;
    int value;
};

struct Template {
    char label;
    unsigned int id;
    std::vector<coord> coords;
};

struct template_entry {
    char label;
    unsigned int id;
    int x;
    int y;
    int z;
    int value;
};

template_entry process_entry(std::string line, int line_num) {
    template_entry e;
    std::stringstream ss(line);                     // create a string stream object to parse the line
    char comma;
    if (!(ss >> e.label >> comma >> e.id >> comma >> e.x >> comma >> e.y >> comma >> e.z >> comma >> e.value)) {
        throw std::runtime_error("Could not parse line " + line_num);
    }
    return e;
}

// load templates from CSV
std::vector<Template> load_templates_from_csv(const std::string& filename) {
    std::ifstream file(filename);                                                   // open the template csv file for reading
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::vector<Template> templates;                                                // create an empty list of templates
    //std::map<std::pair<char, int>, Template> template_map; // to store (class, ID) -> Template



    std::string line;                               // create a line used to buffer a line of text in the input file
    std::getline(file, line);                   // throw away the first line

    unsigned int line_number = 0;
    unsigned int template_number = 0;
    Template t;                                         // temporary template used to store values while reading
    while (std::getline(file, line)) {              // while there are still lines in the file, read a line
        line_number++;                                  // increment the line number

        template_entry e = process_entry(line, line_number);    // process an entry from the CSV file
        if (line_number == 1) {                                     // if this is the first line
            t.label = e.label;                                      // create the first template
            t.id = e.id;
            t.coords.clear();
        }
        else if (t.label != e.label || t.id != e.id) {              // if the current line represents a new template
            templates.push_back(t);                                 // store the previous template
            t.label = e.label;                                      // update the label and ID
            t.id = e.id;
            t.coords.clear();                                       // empty the coordinates
        }
        coord c;                                                    // create and fill the coordinate represented by this line
        c.dx = e.x;
        c.dy = e.y;
        c.dz = e.z;
        c.value = e.value;  //  assign the value field!
        t.coords.push_back(c);                                      // add the coordinate to the template

    }

    // map values to vector
    //for (const auto& pair : template_map) {
    //    templates.push_back(pair.second);
    //}

    return templates;
}



struct Point {
    int x, y, z;
};

// Initialize a zero-volume
tira::volume<int> zero_volume(int x_size, int y_size, int z_size) {
    tira::volume<int> volume(x_size, y_size, z_size);
    for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
            for (int z = 0; z < z_size; z++) {
                volume(x, y, z) = 0;
            }
        }
    }
    return volume;
}

bool is_p_simple(tira::volume<int>& volume, int x, int y, int z) {
    int T26 = 0, T6 = 0;
    std::vector<Point> neighbors_26, neighbors_6 = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}
    };

    // background R T26 (volume == 0)
    std::vector<Point> R_neighbors;
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            for (int dz = -1; dz <= 1; dz++) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                int nx = x + dx, ny = y + dy, nz = z + dz;
                if (nx >= 0 && nx < volume.X() && ny >= 0 && ny < volume.Y() && nz >= 0 && nz < volume.Z()) {
                    if (volume(nx, ny, nz) == 255) {
                        T26++;
                    }
                    else {
                        R_neighbors.push_back({ nx, ny, nz });
                    }
                }
            }
        }
    }

    // Compute T6 (Background Connectivity)
    for (const auto& offset : neighbors_6) {
        int nx = x + offset.x, ny = y + offset.y, nz = z + offset.z;
        if (nx >= 0 && nx < volume.X() && ny >= 0 && ny < volume.Y() && nz >= 0 && nz < volume.Z()) {
            if (volume(nx, ny, nz) == 0) {
                T6++;
            }
        }
    }

    //  adjacency condition
    for (const auto& y : R_neighbors) {
        bool found = false;
        for (const auto& z : R_neighbors) {
            if ((abs(z.x - x) <= 1 && abs(z.y - y.y) <= 1 && abs(z.z - y.z) <= 1) &&
                (abs(z.x - y.x) <= 1 && abs(z.y - y.y) <= 1 && abs(z.z - y.z) <= 1)) {
                found = true;
                break;
            }
        }
        if (!found) return false; // Not P-simple
    }


    // Corrected P-simple condition
    return (T26 > 1 && T6 > 1);
}




// Check if voxel matches any template (same as MS96)
bool matches_template(tira::volume<int>& volume, int x, int y, int z, const Template& tmpl) {
    bool has_at_least_one_1 = false;
    int x_count = 0;

    for (const auto& coord : tmpl.coords) {
        int nx = x + coord.dx, ny = y + coord.dy, nz = z + coord.dz;

        if (nx < 0 || nx >= volume.X() || ny < 0 || ny >= volume.Y() || nz < 0 || nz >= volume.Z()) {
            if (coord.value != -1 && coord.value != 0) return false;
            continue;
        }

        int voxel_value = volume(nx, ny, nz);

        if (coord.value == -1) {
            x_count++;
            if (voxel_value == 255) {
                has_at_least_one_1 = true;
            }
        }
        else if (voxel_value != coord.value) {
            return false;
        }
    }

    if (x_count > 0) return has_at_least_one_1;
    return true;
}

bool is_tail_point(tira::volume<int>& volume, int x, int y, int z) {
    int neighbor_count = 0;
    std::vector<Point> neighbors_26 = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
        {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
        {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
        {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
        {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}
    };

    //  number of 26-neighbors that are in the object (value = 255)
    for (const auto& offset : neighbors_26) {
        int nx = x + offset.x, ny = y + offset.y, nz = z + offset.z;
        if (nx >= 0 && nx < volume.X() && ny >= 0 && ny < volume.Y() && nz >= 0 && nz < volume.Z()) {
            if (volume(nx, ny, nz) == 255) {
                neighbor_count++;
            }
        }
    }

    // Check if it is a tail point
    return (neighbor_count == 1 || neighbor_count == 2); // Line-end or near-line-end
}


// Thinning function P-simple point check
tira::volume<int> thinning_lohou(tira::volume<int> volume, const std::vector<Template>& templates) {
    tira::volume<int> output = volume;
    bool changed;

    do {
        changed = false;
        tira::volume<int> to_delete = zero_volume(volume.X(), volume.Y(), volume.Z());

        for (int y = 0; y < volume.Y(); y++) {
            for (int x = 0; x < volume.X(); x++) {
                for (int z = 0; z < volume.Z(); z++) {
                    if (output(x, y, z) == 255) {
                        if (!is_tail_point(output, x, y, z)) { // **Skip tail points**
                            for (const auto& tmpl : templates) {
                                if (matches_template(output, x, y, z, tmpl) && is_p_simple(output, x, y, z)) {
                                    to_delete(x, y, z) = 255;
                                    changed = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }


        for (int y = 0; y < volume.Y(); y++) {
            for (int x = 0; x < volume.X(); x++) {
                for (int z = 0; z < volume.Z(); z++) {
                    if (to_delete(x, y, z) == 255) {
                        output(x, y, z) = 0;
                    }
                }
            }
        }
    } while (changed);

    return output;
}

// Lohou thinning function that directly works with tira::volume<int>
void lohou(tira::volume<int>& in, tira::volume<int>& out, int x, int y, int z, const std::string& csv_file) {

    // Load templates from the CSV file 
    std::vector<Template> templates = load_templates_from_csv(csv_file);

    out = tira::volume<int>(x, y, z);

    //  Binarize 
    /*for (int zi = 0; zi < z; ++zi) {
        for (int yi = 0; yi < y; ++yi) {
            for (int xi = 0; xi < x; ++xi) {
                if (in(xi, yi, zi) != 0)
                    in(xi, yi, zi) = 1;
                else
                    in(xi, yi, zi) = 0;
            }
        }
    }*/

    // Lohou's 3D thinning 
    thinning_lohou(in, templates);

    // Copy the result to the output volume
    for (int zi = 0; zi < z; ++zi) {
        for (int yi = 0; yi < y; ++yi) {
            for (int xi = 0; xi < x; ++xi) {
                out(xi, yi, zi) = in(xi, yi, zi);
            }
        }
    }
}