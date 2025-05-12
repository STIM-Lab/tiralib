//void rsf3d(float* phi, float* image, float* Force, float* Regularization, float* Smoothness, float* Band, unsigned long sx, unsigned long sy, unsigned long sz, float sigma_r, float sigma_F, float epsilon);




#include <tira/volume.h>
#include <vector>
#include <array>
#include <set>
#include <queue>
#include <tuple>
#include <iostream>
#include <algorithm>
#include <numeric>


using Volume = tira::volume<int>;
using Point = std::array<int, 3>;
using PointList = std::vector<Point>;


void prepare_data(Volume& volume) {
    for (int y = 0; y < volume.Y(); ++y)
        for (int x = 0; x < volume.X(); ++x)
            for (int z = 0; z < volume.Z(); ++z)
                if (volume(x, y, z) != 0)
                    volume(x, y, z) = 1;
}


// to check out of boundary 
int get_pixel(Volume& vol, int x, int y, int z) {
    if (x >= 0 && x < vol.X() && y >= 0 && y < vol.Y() && z >= 0 && z < vol.Z()) {
        return vol(x, y, z);
    }
    return 0;
}


/* -----------------------------------------------------------------------*/
    /**
     * Get neighborhood of a pixel in a 3D image (0 border conditions)
     * Retrieves the full 3x3x3 voxel neighborhood of a given voxel (including center)
     * @param image 3D image
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @return corresponding 27-pixels neighborhood (0 if out of image)
     */

std::array<int, 27> get_neighborhood(Volume& vol, int x, int y, int z) {
    std::array<int, 27> neighborhood;
    neighborhood[0] = get_pixel(vol, x - 1, y - 1, z - 1);
    neighborhood[1] = get_pixel(vol, x, y - 1, z - 1);
    neighborhood[2] = get_pixel(vol, x + 1, y - 1, z - 1);

    neighborhood[3] = get_pixel(vol, x - 1, y, z - 1);
    neighborhood[4] = get_pixel(vol, x, y, z - 1);
    neighborhood[5] = get_pixel(vol, x + 1, y, z - 1);

    neighborhood[6] = get_pixel(vol, x - 1, y + 1, z - 1);
    neighborhood[7] = get_pixel(vol, x, y + 1, z - 1);
    neighborhood[8] = get_pixel(vol, x + 1, y + 1, z - 1);

    neighborhood[9] = get_pixel(vol, x - 1, y - 1, z);
    neighborhood[10] = get_pixel(vol, x, y - 1, z);
    neighborhood[11] = get_pixel(vol, x + 1, y - 1, z);

    neighborhood[12] = get_pixel(vol, x - 1, y, z);
    neighborhood[13] = get_pixel(vol, x, y, z);
    neighborhood[14] = get_pixel(vol, x + 1, y, z);

    neighborhood[15] = get_pixel(vol, x - 1, y + 1, z);
    neighborhood[16] = get_pixel(vol, x, y + 1, z);
    neighborhood[17] = get_pixel(vol, x + 1, y + 1, z);

    neighborhood[18] = get_pixel(vol, x - 1, y - 1, z + 1);
    neighborhood[19] = get_pixel(vol, x, y - 1, z + 1);
    neighborhood[20] = get_pixel(vol, x + 1, y - 1, z + 1);

    neighborhood[21] = get_pixel(vol, x - 1, y, z + 1);
    neighborhood[22] = get_pixel(vol, x, y, z + 1);
    neighborhood[23] = get_pixel(vol, x + 1, y, z + 1);

    neighborhood[24] = get_pixel(vol, x - 1, y + 1, z + 1);
    neighborhood[25] = get_pixel(vol, x, y + 1, z + 1);
    neighborhood[26] = get_pixel(vol, x + 1, y + 1, z + 1);
    return neighborhood;
}

//it only connects to one other voxel in the neighborhood
//* Check if a point in the given stack is at the end of an arc
//true if the point has exactly one neighbor
bool is_endpoint(Volume& vol, int x, int y, int z) {
    auto neighbor = get_neighborhood(vol, x, y, z);
    int count = -1;
    for (int i = 0; i < 27; ++i) {
        if (neighbor[i] == 1) count++;
    }
    return count == 1;
}

/*
* Precomputes a lookup table of size 256 (for all possible 8-bit combinations of a 2x2x2 voxel cube).
Each entry encodes the contribution to the Euler characteristic.
/**
     * Fill Euler LUT
     *
     * @param LUT Euler LUT

*/
std::array<int, 256> fill_euler_LUT() {
    std::array<int, 256> LUT = { 0 };

    LUT[1] = 1;  LUT[3] = -1; LUT[5] = -1; LUT[7] = 1; LUT[9] = -3; LUT[11] = -1;
    LUT[13] = -1; LUT[15] = 1; LUT[17] = -1; LUT[19] = 1; LUT[21] = 1; LUT[23] = -1;
    LUT[25] = 3; LUT[27] = 1; LUT[29] = 1; LUT[31] = -1; LUT[33] = -3; LUT[35] = -1;
    LUT[37] = 3; LUT[39] = 1; LUT[41] = 1; LUT[43] = -1; LUT[45] = 3; LUT[47] = 1;
    LUT[49] = -1; LUT[51] = 1; LUT[53] = 1; LUT[55] = -1; LUT[57] = 3; LUT[59] = 1;
    LUT[61] = 1; LUT[63] = -1; LUT[65] = -3; LUT[67] = 3; LUT[69] = -1; LUT[71] = 1;
    LUT[73] = 1; LUT[75] = 3; LUT[77] = -1; LUT[79] = 1; LUT[81] = -1; LUT[83] = 1;
    LUT[85] = 1; LUT[87] = -1; LUT[89] = 3; LUT[91] = 1; LUT[93] = 1; LUT[95] = -1;
    LUT[97] = 1; LUT[99] = 3; LUT[101] = 3; LUT[103] = 1; LUT[105] = 5; LUT[107] = 3;
    LUT[109] = 3; LUT[111] = 1; LUT[113] = -1; LUT[115] = 1; LUT[117] = 1; LUT[119] = -1;
    LUT[121] = 3; LUT[123] = 1; LUT[125] = 1; LUT[127] = -1; LUT[129] = -7; LUT[131] = -1;
    LUT[133] = -1; LUT[135] = 1; LUT[137] = -3; LUT[139] = -1; LUT[141] = -1; LUT[143] = 1;
    LUT[145] = -1; LUT[147] = 1; LUT[149] = 1; LUT[151] = -1; LUT[153] = 3; LUT[155] = 1;
    LUT[157] = 1; LUT[159] = -1; LUT[161] = -3; LUT[163] = -1; LUT[165] = 3; LUT[167] = 1;
    LUT[169] = 1; LUT[171] = -1; LUT[173] = 3; LUT[175] = 1; LUT[177] = -1; LUT[179] = 1;
    LUT[181] = 1; LUT[183] = -1; LUT[185] = 3; LUT[187] = 1; LUT[189] = 1; LUT[191] = -1;
    LUT[193] = -3; LUT[195] = 3; LUT[197] = -1; LUT[199] = 1; LUT[201] = 1; LUT[203] = 3;
    LUT[205] = -1; LUT[207] = 1; LUT[209] = -1; LUT[211] = 1; LUT[213] = 1; LUT[215] = -1;
    LUT[217] = 3; LUT[219] = 1; LUT[221] = 1; LUT[223] = -1; LUT[225] = 1; LUT[227] = 3;
    LUT[229] = 3; LUT[231] = 1; LUT[233] = 5; LUT[235] = 3; LUT[237] = 3; LUT[239] = 1;
    LUT[241] = -1; LUT[243] = 1; LUT[245] = 1; LUT[247] = -1; LUT[249] = 3; LUT[251] = 1;
    LUT[253] = 1; LUT[255] = -1;

    return LUT;
}

//Directional accessors to test border types based on whether adjacent voxels in that direction are background 

int N(Volume& vol, int x, int y, int z) { return get_pixel(vol, x, y - 1, z); }
int S(Volume& vol, int x, int y, int z) { return get_pixel(vol, x, y + 1, z); }
int E(Volume& vol, int x, int y, int z) { return get_pixel(vol, x + 1, y, z); }
int W(Volume& vol, int x, int y, int z) { return get_pixel(vol, x - 1, y, z); }
int U(Volume& vol, int x, int y, int z) { return get_pixel(vol, x, y, z + 1); }
int B(Volume& vol, int x, int y, int z) { return get_pixel(vol, x, y, z - 1); }





int get_pixel_nocheck(Volume& vol, int x, int y, int z) {
    return vol(x, y, z);
}

void set_pixel(Volume& vol, int x, int y, int z, int value) {
    if (x >= 0 && x < vol.X() && y >= 0 && y < vol.Y() && z >= 0 && z < vol.Z()) {
        vol(x, y, z) = value;
    }
}




//Creates a LUT that maps integers 0–255 to their number of 1-bits (i.e., how many neighbors are on). 
// This isn’t always used, but it helps quickly evaluate connectivity degree.
//Fill number of points in octant LUT
void fill_num_of_points_LUT(std::array<int, 256>& LUT) {
    for (int i = 0; i < 256; ++i) {
        int count = 0;
        int val = i;
        while (val) {
            count += val & 1;
            val >>= 1;
        }
        LUT[i] = count;
    }
}


//Each function computes a bitmask index from a specific octant of the 3x3x3 neighborhood. 
// These 7-voxel octants are used to query the Euler LUT to ensure topological invariance. 
// The bits encode which of the 7 voxels in the octant are present (1) or absent (0).
/**
     * Check if a point is Euler invariant
     *
     * @param neighbors neighbor pixels of the point
     * @param LUT Euler LUT
     * @return true or false if the point is Euler invariant or not
     */

uint8_t index_octant_NEB(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[2])   v |= 128;
    if (n[1])   v |= 64;
    if (n[11])  v |= 32;
    if (n[10])  v |= 16;
    if (n[5])   v |= 8;
    if (n[4])   v |= 4;
    if (n[14])  v |= 2;
    return v;
}

uint8_t index_octant_NWB(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[0])   v |= 128;
    if (n[9])   v |= 64;
    if (n[3])   v |= 32;
    if (n[12])  v |= 16;
    if (n[1])   v |= 8;
    if (n[10])  v |= 4;
    if (n[4])   v |= 2;
    return v;
}

uint8_t index_octant_SEB(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[8])   v |= 128;
    if (n[7])   v |= 64;
    if (n[17])  v |= 32;
    if (n[16])  v |= 16;
    if (n[5])   v |= 8;
    if (n[4])   v |= 4;
    if (n[14])  v |= 2;
    return v;
}

uint8_t index_octant_SWB(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[6])   v |= 128;
    if (n[15])  v |= 64;
    if (n[7])   v |= 32;
    if (n[16])  v |= 16;
    if (n[3])   v |= 8;
    if (n[12])  v |= 4;
    if (n[4])   v |= 2;
    return v;
}

uint8_t index_octant_NEU(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[20])  v |= 128;
    if (n[23])  v |= 64;
    if (n[19])  v |= 32;
    if (n[22])  v |= 16;
    if (n[11])  v |= 8;
    if (n[14])  v |= 4;
    if (n[10])  v |= 2;
    return v;
}

uint8_t index_octant_NWU(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[18]) v |= 128;
    if (n[21]) v |= 64;
    if (n[9])  v |= 32;
    if (n[12]) v |= 16;
    if (n[19]) v |= 8;
    if (n[22]) v |= 4;
    if (n[10]) v |= 2;
    return v;
}

uint8_t index_octant_SEU(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[26]) v |= 128;
    if (n[23]) v |= 64;
    if (n[17]) v |= 32;
    if (n[14]) v |= 16;
    if (n[25]) v |= 8;
    if (n[22]) v |= 4;
    if (n[16]) v |= 2;
    return v;
}

uint8_t index_octant_SWU(const std::array<uint8_t, 27>& n) {
    uint8_t v = 1;
    if (n[24]) v |= 128;
    if (n[25]) v |= 64;
    if (n[15]) v |= 32;
    if (n[16]) v |= 16;
    if (n[21]) v |= 8;
    if (n[22]) v |= 4;
    if (n[12]) v |= 2;
    return v;
}

//Checks all 8 octants using their corresponding index_octant_* functions and sums the Euler LUT values. 
// If the sum is 0, the deletion preserves Euler characteristic.

bool is_euler_invariant(const std::array<uint8_t, 27>& neighbors, const std::array<int, 256>& LUT) {
    int eulerChar = 0;
    eulerChar += LUT[index_octant_SWU(neighbors)];
    eulerChar += LUT[index_octant_SEU(neighbors)];
    eulerChar += LUT[index_octant_NWU(neighbors)];
    eulerChar += LUT[index_octant_NEU(neighbors)];
    eulerChar += LUT[index_octant_SWB(neighbors)];
    eulerChar += LUT[index_octant_SEB(neighbors)];
    eulerChar += LUT[index_octant_NWB(neighbors)];
    eulerChar += LUT[index_octant_NEB(neighbors)];
    return eulerChar == 0;
}


// recursive labeling function that marks all connected voxels in a given octant with the same label. 
// It simulates connected component labeling within each octant.
/* -----------------------------------------------------------------------*/
    /**
     * This is a recursive method that calculates the number of connected
     * components in the 3D neighborhood after the center pixel would
     * have been removed.
     *
     * @param octant
     * @param label
     * @param cube
     */


void octree_labeling(int octant, int label, std::array<int, 26>& cube) {
    if (octant == 1) {
        if (cube[0] == 1) cube[0] = label;
        if (cube[1] == 1) { cube[1] = label; octree_labeling(2, label, cube); }
        if (cube[3] == 1) { cube[3] = label; octree_labeling(3, label, cube); }
        if (cube[4] == 1) {
            cube[4] = label;
            octree_labeling(2, label, cube);
            octree_labeling(3, label, cube);
            octree_labeling(4, label, cube);
        }
        if (cube[9] == 1) { cube[9] = label; octree_labeling(5, label, cube); }
        if (cube[10] == 1) {
            cube[10] = label;
            octree_labeling(2, label, cube);
            octree_labeling(5, label, cube);
            octree_labeling(6, label, cube);
        }
        if (cube[12] == 1) {
            cube[12] = label;
            octree_labeling(3, label, cube);
            octree_labeling(5, label, cube);
            octree_labeling(7, label, cube);
        }
    }
    if (octant == 2) {
        if (cube[1] == 1) { cube[1] = label; octree_labeling(1, label, cube); }
        if (cube[4] == 1) {
            cube[4] = label;
            octree_labeling(1, label, cube);
            octree_labeling(3, label, cube);
            octree_labeling(4, label, cube);
        }
        if (cube[10] == 1) {
            cube[10] = label;
            octree_labeling(1, label, cube);
            octree_labeling(5, label, cube);
            octree_labeling(6, label, cube);
        }
        if (cube[2] == 1) cube[2] = label;
        if (cube[5] == 1) { cube[5] = label; octree_labeling(4, label, cube); }
        if (cube[11] == 1) { cube[11] = label; octree_labeling(6, label, cube); }
        if (cube[13] == 1) {
            cube[13] = label;
            octree_labeling(4, label, cube);
            octree_labeling(6, label, cube);
            octree_labeling(8, label, cube);
        }
    }
    if (octant == 3) {
        if (cube[3] == 1) { cube[3] = label; octree_labeling(1, label, cube); }
        if (cube[4] == 1) {
            cube[4] = label;
            octree_labeling(1, label, cube);
            octree_labeling(2, label, cube);
            octree_labeling(4, label, cube);
        }
        if (cube[12] == 1) {
            cube[12] = label;
            octree_labeling(1, label, cube);
            octree_labeling(5, label, cube);
            octree_labeling(7, label, cube);
        }
        if (cube[6] == 1) cube[6] = label;
        if (cube[7] == 1) { cube[7] = label; octree_labeling(4, label, cube); }
        if (cube[14] == 1) { cube[14] = label; octree_labeling(7, label, cube); }
        if (cube[15] == 1) {
            cube[15] = label;
            octree_labeling(4, label, cube);
            octree_labeling(7, label, cube);
            octree_labeling(8, label, cube);
        }
    }
    if (octant == 4) {
        if (cube[4] == 1) {
            cube[4] = label;
            octree_labeling(1, label, cube);
            octree_labeling(2, label, cube);
            octree_labeling(3, label, cube);
        }
        if (cube[5] == 1) { cube[5] = label; octree_labeling(2, label, cube); }
        if (cube[13] == 1) {
            cube[13] = label;
            octree_labeling(2, label, cube);
            octree_labeling(6, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[7] == 1) { cube[7] = label; octree_labeling(3, label, cube); }
        if (cube[15] == 1) {
            cube[15] = label;
            octree_labeling(3, label, cube);
            octree_labeling(7, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[8] == 1) cube[8] = label;
        if (cube[16] == 1) { cube[16] = label; octree_labeling(8, label, cube); }
    }
    if (octant == 5) {
        if (cube[9] == 1) { cube[9] = label; octree_labeling(1, label, cube); }
        if (cube[10] == 1) {
            cube[10] = label;
            octree_labeling(1, label, cube);
            octree_labeling(2, label, cube);
            octree_labeling(6, label, cube);
        }
        if (cube[12] == 1) {
            cube[12] = label;
            octree_labeling(1, label, cube);
            octree_labeling(3, label, cube);
            octree_labeling(7, label, cube);
        }
        if (cube[17] == 1) cube[17] = label;
        if (cube[18] == 1) { cube[18] = label; octree_labeling(6, label, cube); }
        if (cube[20] == 1) { cube[20] = label; octree_labeling(7, label, cube); }
        if (cube[21] == 1) {
            cube[21] = label;
            octree_labeling(6, label, cube);
            octree_labeling(7, label, cube);
            octree_labeling(8, label, cube);
        }
    }
    if (octant == 6) {
        if (cube[10] == 1) {
            cube[10] = label;
            octree_labeling(1, label, cube);
            octree_labeling(2, label, cube);
            octree_labeling(5, label, cube);
        }
        if (cube[11] == 1) { cube[11] = label; octree_labeling(2, label, cube); }
        if (cube[13] == 1) {
            cube[13] = label;
            octree_labeling(2, label, cube);
            octree_labeling(4, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[18] == 1) { cube[18] = label; octree_labeling(5, label, cube); }
        if (cube[21] == 1) {
            cube[21] = label;
            octree_labeling(5, label, cube);
            octree_labeling(7, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[19] == 1) cube[19] = label;
        if (cube[22] == 1) { cube[22] = label; octree_labeling(8, label, cube); }
    }
    if (octant == 7) {
        if (cube[12] == 1) {
            cube[12] = label;
            octree_labeling(1, label, cube);
            octree_labeling(3, label, cube);
            octree_labeling(5, label, cube);
        }
        if (cube[14] == 1) { cube[14] = label; octree_labeling(3, label, cube); }
        if (cube[15] == 1) {
            cube[15] = label;
            octree_labeling(3, label, cube);
            octree_labeling(4, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[20] == 1) { cube[20] = label; octree_labeling(5, label, cube); }
        if (cube[21] == 1) {
            cube[21] = label;
            octree_labeling(5, label, cube);
            octree_labeling(6, label, cube);
            octree_labeling(8, label, cube);
        }
        if (cube[23] == 1) cube[23] = label;
        if (cube[24] == 1) { cube[24] = label; octree_labeling(8, label, cube); }
    }
    if (octant == 8) {
        if (cube[13] == 1) {
            cube[13] = label;
            octree_labeling(2, label, cube);
            octree_labeling(4, label, cube);
            octree_labeling(6, label, cube);
        }
        if (cube[15] == 1) {
            cube[15] = label;
            octree_labeling(3, label, cube);
            octree_labeling(4, label, cube);
            octree_labeling(7, label, cube);
        }
        if (cube[16] == 1) { cube[16] = label; octree_labeling(4, label, cube); }
        if (cube[21] == 1) {
            cube[21] = label;
            octree_labeling(5, label, cube);
            octree_labeling(6, label, cube);
            octree_labeling(7, label, cube);
        }
        if (cube[22] == 1) { cube[22] = label; octree_labeling(6, label, cube); }
        if (cube[24] == 1) { cube[24] = label; octree_labeling(7, label, cube); }
        if (cube[25] == 1) cube[25] = label;
    }
}


//Determines whether a voxel is “simple” which means., its removal doesn’t break connectivity. 
// It copies the 26-neighborhood (skipping center), performs labeling via octree_labeling, 
// and checks if only one connected component exists. 
// If multiple are found, it’s not simple.
/**
     * Check if current point is a Simple Point.
     * This method is named 'N(v)_labeling' in [Lee94].
     * Outputs the number of connected objects in a neighborhood of a point
     * after this point would have been removed.
     *
     * @param neighbors neighbor pixels of the point
     * @return true or false if the point is simple or not
     */
bool is_simple_point(const std::array<uint8_t, 27>& neighbors) {
    std::array<int, 26> cube;

    // Copy first 13 elements as-is
    for (int i = 0; i < 13; ++i)
        cube[i] = neighbors[i];

    // Skip index 13 (center), then continue from 14 to 26 into cube[13 to 25]
    for (int i = 14; i < 27; ++i)
        cube[i - 1] = neighbors[i];

    int label = 2;

    for (int i = 0; i < 26; ++i) {
        if (cube[i] == 1) {
            switch (i) {
            case 0: case 1: case 3: case 4: case 9: case 10: case 12:
                octree_labeling(1, label, cube);
                break;
            case 2: case 5: case 11: case 13:
                octree_labeling(2, label, cube);
                break;
            case 6: case 7: case 14: case 15:
                octree_labeling(3, label, cube);
                break;
            case 8: case 16:
                octree_labeling(4, label, cube);
                break;
            case 17: case 18: case 20: case 21:
                octree_labeling(5, label, cube);
                break;
            case 19: case 22:
                octree_labeling(6, label, cube);
                break;
            case 23: case 24:
                octree_labeling(7, label, cube);
                break;
            case 25:
                octree_labeling(8, label, cube);
                break;
            }
            label++;
            if ((label - 2) >= 2) {
                return false;
            }
        }
    }

    return true;
}

/*
* Main function that performs thinning by iterating through each of the 6 borders. For each voxel, it checks:

Is it a border voxel (based on direction)?

Is it not an endpoint?

Is it Euler-invariant?

Is it simple?

If all conditions are met, the voxel is queued for deletion. After each pass,
a second pass confirms deletability to prevent conflicts, then deletion proceeds.
This loop continues until no voxels are deleted in 6 successive directional passes.


*/


void computeThinImage(Volume& volume) {
    int width = volume.X();
    int height = volume.Y();
    int depth = volume.Z();

    std::array<int, 256> eulerLUT = fill_euler_LUT();
    std::array<int, 256> pointsLUT;
    fill_num_of_points_LUT(pointsLUT);

    std::vector<Point> simpleBorderPoints;
    int iterations = 0;
    int unchangedBorders = 0;

    while (unchangedBorders < 6) {
        unchangedBorders = 0;
        iterations++;

        for (int currentBorder = 1; currentBorder <= 6; currentBorder++) {
            bool noChange = true;

            for (int z = 0; z < depth; z++) {
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        if (get_pixel_nocheck(volume, x, y, z) != 1)
                            continue;

                        bool isBorderPoint = false;

                        if (currentBorder == 1 && N(volume, x, y, z) <= 0) isBorderPoint = true;
                        if (currentBorder == 2 && S(volume, x, y, z) <= 0) isBorderPoint = true;
                        if (currentBorder == 3 && E(volume, x, y, z) <= 0) isBorderPoint = true;
                        if (currentBorder == 4 && W(volume, x, y, z) <= 0) isBorderPoint = true;
                        if (currentBorder == 5 && U(volume, x, y, z) <= 0) isBorderPoint = true;
                        if (currentBorder == 6 && B(volume, x, y, z) <= 0) isBorderPoint = true;

                        if (!isBorderPoint)
                            continue;

                        if (is_endpoint(volume, x, y, z))
                            continue;

                        std::array<int, 27> neighborhood_int = get_neighborhood(volume, x, y, z);
                        std::array<uint8_t, 27> neighborhood;
                        for (int i = 0; i < 27; ++i) neighborhood[i] = static_cast<uint8_t>(neighborhood_int[i]);

                        if (!is_euler_invariant(neighborhood, eulerLUT))
                            continue;

                        if (!is_simple_point(neighborhood))
                            continue;


                        simpleBorderPoints.push_back({ x, y, z });
                    }
                }
            }

            for (const auto& index : simpleBorderPoints) {
                std::array<int, 27> neighbors_int = get_neighborhood(volume, index[0], index[1], index[2]);
                std::array<uint8_t, 27> neighbors;
                for (int i = 0; i < 27; ++i) neighbors[i] = static_cast<uint8_t>(neighbors_int[i]);

                if (is_simple_point(neighbors)) {
                    set_pixel(volume, index[0], index[1], index[2], 0);
                    noChange = false;
                }

            }

            if (noChange)
                unchangedBorders++;

            simpleBorderPoints.clear();
        }
    }
}

// Lee thinning function that directly works with tira::volume<int>
void lee(tira::volume<int>& in, tira::volume<int>& out, int x, int y, int z) {
    
    out = tira::volume<int>(x, y, z);

    //  Binarize 
    for (int zi = 0; zi < z; ++zi) {
        for (int yi = 0; yi < y; ++yi) {
            for (int xi = 0; xi < x; ++xi) {
                if (in(xi, yi, zi) != 0)
                    in(xi, yi, zi) = 1;
                else
                    in(xi, yi, zi) = 0;
            }
        }
    }

    // Lee's 3D thinning 
    computeThinImage(in);

    // Copy the result to the output volume
    for (int zi = 0; zi < z; ++zi) {
        for (int yi = 0; yi < y; ++yi) {
            for (int xi = 0; xi < x; ++xi) {
                out(xi, yi, zi) = in(xi, yi, zi);
            }
        }
    }
}

