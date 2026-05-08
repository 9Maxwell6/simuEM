#pragma once
#include <bit>

#include "utils/logger.h"
#include "utils/util_hash.h"


namespace simu {


using namespace util;


struct Hash_ID
{
    size_t id;
    bool exist;
};


struct Hash_Node
{
private:
    size_t node_id = 0;
    std::vector<Block_1> table_1;
    // mask for allocating block entry
    size_t mask_1;
    // increment on each insert, for generating new id and rehash
    size_t count_1 = 0;
    void rehash_1(size_t new_size);

public:
    Hash_Node(size_t initial_size_1=0);

    size_t get_id(size_t p0, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p_dof);

    const size_t get_node_dof_count() const { return count_1; }

    void clear_table();
};


struct Hash_Edge
{
private:
    size_t edge_id = 0;
    std::vector<Block_2> table_2;
    size_t mask_2;
    size_t count_2 = 0;
    void rehash_2(size_t new_size);

public:
    Hash_Edge(size_t initial_size_2=0);

    size_t get_id(size_t p0, size_t p1, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p1, size_t p_dof);

    const size_t get_edge_dof_count() const { return count_2; }

    void clear_table();
};



struct Hash_Face
{
private:
    size_t face_id = 0;
    std::vector<Block_3> table_3;   // for 3 node face (e.g. face of tetrahedron).
    size_t mask_3;
    size_t count_3 = 0;
    void rehash_3(size_t new_size);

    std::vector<Block_4> table_4;   // for 4 node face (e.g. face of hexahedron).
    size_t mask_4;
    size_t count_4 = 0;
    void rehash_4(size_t new_size);

public:
    Hash_Face(size_t initial_size_3=0, size_t initial_size_4=0);

    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p1, size_t p2, size_t p_dof);

    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof);

    const size_t get_face_dof_count() const { return count_3+count_4; }

    void clear_table();
};

struct Hash_Cell
{
private:
    size_t cell_id = 0;
    std::vector<Block_3> table_3;   // for 3 node cell (triangle).
    size_t mask_3;
    size_t count_3 = 0;
    void rehash_3(size_t new_size);

    std::vector<Block_4> table_4;   // for 4 node cell (tetrahedron).
    size_t mask_4;
    size_t count_4 = 0;
    void rehash_4(size_t new_size);


public:
    Hash_Cell(size_t initial_size_3=0, size_t initial_size_4=0);

    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p1, size_t p2, size_t p_dof);

    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof);
    Hash_ID exist(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof);

    const size_t get_cell_dof_count() const { return count_3+count_4; }

    void clear_table();
};



}