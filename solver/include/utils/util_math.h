#pragma once

#include "entity/mesh/node.h"

#include <vector>
#include <algorithm>




namespace util 
{

/**
 * @brief Return a pair of two ids sorted in ascending order.
 * 
 * @param a First id.
 * @param b Second id.
 * @return A vector of size 2 with the smaller id first and the larger id second.
 */
inline std::vector<size_t> sort_a_b(size_t a, size_t b) {
    if (a < b) return {a, b};
    return {b, a};
}




/**
 * @brief Order three triangle node ids so that the face normal points outward,
 *        away from the opposite vertex.
 * 
 * The outward direction is determined by computing the cross product normal
 * of the triangle and checking its orientation against the opposite vertex.
 * If the normal points inward (toward the opposite vertex), the winding
 * order of b and c is swapped to flip the normal.
 * 
 * @param a  node id of the triangle.
 * @param b  node id of the triangle.
 * @param c  node id of the triangle.
 * @param pa coordinate of node a.
 * @param pb coordinate of node b.
 * @param pc coordinate of node c.
 * @param opposite coordinate of the vertex opposite to this face.
 * @return A vector of three node ids with outward-facing winding order.
 */
inline std::vector<size_t> sort_outward_triangle(size_t a, size_t b, size_t c,
                                                 const simu::Node& pa,
                                                 const simu::Node& pb,
                                                 const simu::Node& pc,
                                                 const simu::Node& opposite)
{
    // Edge vectors
    double e1x = pb.x - pa.x, e1y = pb.y - pa.y, e1z = pb.z - pa.z;
    double e2x = pc.x - pa.x, e2y = pc.y - pa.y, e2z = pc.z - pa.z;

    // Cross product (face normal)
    double nx = e1y*e2z - e1z*e2y;
    double ny = e1z*e2x - e1x*e2z;
    double nz = e1x*e2y - e1y*e2x;

    // Dot normal with vector from opposite vertex to pa
    double dot =  nx*(pa.x - opposite.x)
                + ny*(pa.y - opposite.y)
                + nz*(pa.z - opposite.z);

    if (dot < 0) std::swap(b, c);
    
    return {a, b, c};
}


/**
 * @brief Check if none of the elements in vector a exist in vector b.
 * 
 * @param a Vector of ids to check.
 * @param b Vector of ids to check against.
 * @return true if no element in a is found in b, false otherwise.
 */
inline bool a_not_in_b(const std::vector<size_t>& a, const std::vector<size_t>& b)
{
    for (auto id : a)
    {
        if (std::find(b.begin(), b.end(), id) != b.end()) {
            return false;
        }
    }
    return true;
}


}
