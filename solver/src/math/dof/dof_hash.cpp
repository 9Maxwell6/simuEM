#include "math/dof/dof_hash.h"


using namespace simu;




Hash_Node::Hash_Node(size_t initial_size_1) : table_1(initial_size_1), mask_1(initial_size_1 - 1)
{
    if((initial_size_1 & (initial_size_1 - 1)) != 0){
        initial_size_1 = std::bit_ceil(initial_size_1);
        table_1.resize(initial_size_1);
        mask_1 = initial_size_1 - 1;
    }
}

Hash_Edge::Hash_Edge(size_t initial_size_2) : table_2(initial_size_2), mask_2(initial_size_2 - 1)
{
    if((initial_size_2 & (initial_size_2 - 1)) != 0){
        initial_size_2 = std::bit_ceil(initial_size_2);
        table_2.resize(initial_size_2);
        mask_2 = initial_size_2 - 1;
    }
}

Hash_Face::Hash_Face(size_t initial_size_3, size_t initial_size_4) : table_3(initial_size_3), mask_3(initial_size_3 - 1),
                                                                     table_4(initial_size_4), mask_4(initial_size_4 - 1)
{
    if((initial_size_3 & (initial_size_3 - 1)) != 0){
        initial_size_3 = std::bit_ceil(initial_size_3);
        table_3.resize(initial_size_3);
        mask_3 = initial_size_3 - 1;
    }

    if((initial_size_4 & (initial_size_4 - 1)) != 0){
        initial_size_4 = std::bit_ceil(initial_size_4);
        table_4.resize(initial_size_4);
        mask_4 = initial_size_4 - 1;
    }
}

Hash_Cell::Hash_Cell(size_t initial_size_3, size_t initial_size_4) : table_3(initial_size_3), mask_3(initial_size_3 - 1),
                                                                     table_4(initial_size_4), mask_4(initial_size_4 - 1)
{
    if((initial_size_3 & (initial_size_3 - 1)) != 0){
        initial_size_3 = std::bit_ceil(initial_size_3);;
        table_3.resize(initial_size_3);
        mask_3 = initial_size_3 - 1;
    }

    if((initial_size_4 & (initial_size_4 - 1)) != 0){
        initial_size_4 = std::bit_ceil(initial_size_4);;
        table_4.resize(initial_size_4);
        mask_4 = initial_size_4 - 1;
    }
}










void Hash_Node::rehash_1(size_t new_size)
{
    std::vector<Block_1> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_1& block : table_1)
        for (Vertex_1& e : block.entry_1) 
            new_table[hash_1(e.v[0], e.dof) & new_mask].entry_1.push_back(e);
    table_1 = std::move(new_table);
    mask_1  = new_mask;
    Logger::info("Hash_Node::rehash_1 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}

void Hash_Edge::rehash_2(size_t new_size)
{
    std::vector<Block_2> new_table(new_size);
    size_t new_mask = new_size - 1;

    for (Block_2& block : table_2)
        for (Vertex_2& e : block.entry_2)
            new_table[hash_2(e.v[0], e.v[1], e.dof) & new_mask].entry_2.push_back(e);
    table_2 = std::move(new_table);
    mask_2  = new_mask;
    Logger::info("Hash_Edge::rehash_2 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}

void Hash_Face::rehash_3(size_t new_size)
{
    std::vector<Block_3> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_3& block : table_3)
        for (Vertex_3& e : block.entry_3)
            new_table[hash_3(e.v[0], e.v[1], e.v[2], e.dof) & new_mask].entry_3.push_back(e);
    table_3 = std::move(new_table);
    mask_3  = new_mask;
    Logger::info("Hash_Face::rehash_3 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}

void Hash_Cell::rehash_3(size_t new_size)
{
    std::vector<Block_3> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_3& block : table_3)
        for (Vertex_3& e : block.entry_3)
            new_table[hash_3(e.v[0], e.v[1], e.v[2], e.dof) & new_mask].entry_3.push_back(e);
    table_3 = std::move(new_table);
    mask_3  = new_mask;
    Logger::info("Hash_Cell::rehash_3 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}

void Hash_Face::rehash_4(size_t new_size)
{
    std::vector<Block_4> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_4& block : table_4)
        for (Vertex_4& e : block.entry_4)
            new_table[hash_4(e.v[0], e.v[1], e.v[2], e.v[3], e.dof) & new_mask].entry_4.push_back(e);
    table_4 = std::move(new_table);
    mask_4  = new_mask;
    Logger::info("Hash_Face::rehash_4 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}

void Hash_Cell::rehash_4(size_t new_size)
{
    std::vector<Block_4> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_4& block : table_4)
        for (Vertex_4& e : block.entry_4)
            new_table[hash_4(e.v[0], e.v[1], e.v[2], e.v[3], e.dof) & new_mask].entry_4.push_back(e);
    table_4 = std::move(new_table);
    mask_4  = new_mask;
    Logger::info("Hash_Cell::rehash_4 - trigger rehash: new size is: "+std::to_string(new_size)+".");
}










size_t Hash_Node::get_id(size_t p0, size_t p_dof)
{
    if (count_1 > 4*table_1.size()) rehash_1(table_1.size() * 2);

    size_t slot = hash_1(p0, p_dof) & mask_1;
    std::vector<Vertex_1>& block_1 = table_1[slot].entry_1;
    for (Vertex_1& e : block_1)
        if (e.v[0] == p0 && e.dof==p_dof)
            return e.id;

    size_t new_id = node_id++;
    block_1.push_back({{p0}, p_dof, new_id});
    count_1++;
    return new_id;
}

size_t Hash_Edge::get_id(size_t p0, size_t p1, size_t p_dof) 
{
    if (count_2 > 4*table_2.size()) rehash_2(table_2.size() * 2);

    if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
    size_t slot = hash_2(p0, p1, p_dof) & mask_2;      
    std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block

    for (Vertex_2& e : block_2)                               // linear scan
        if (e.v[0] == p0 && e.v[1] == p1  && e.dof==p_dof)
            return e.id;

    size_t new_id = edge_id++;                                // global id in table
    block_2.push_back({{p0, p1}, p_dof, new_id});
    count_2++;
    return new_id;
}

size_t Hash_Face::get_id(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    if (count_3 > 4*table_3.size()) rehash_3(table_3.size() * 2);

    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;

    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return e.id;

    size_t new_id = face_id++;  
    block_3.push_back({{v[0], v[1], v[2]}, p_dof, new_id});
    count_3++;
    return new_id;
}

size_t Hash_Cell::get_id(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    if (count_3 > 4*table_3.size()) rehash_3(table_3.size() * 2);

    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;

    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return e.id;

    size_t new_id = cell_id++;  
    block_3.push_back({{v[0], v[1], v[2]}, p_dof, new_id});
    count_3++;
    return new_id;
}

size_t Hash_Face::get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    if (count_4 > 4*table_4.size()) rehash_4(table_4.size() * 2);

    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;

    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return e.id;

    size_t new_id = face_id++;
    block_4.push_back({{v[0], v[1], v[2], v[3]}, p_dof, new_id});
    count_4++;
    return new_id;
}

size_t Hash_Cell::get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    if (count_4 > 4*table_4.size()) rehash_4(table_4.size() * 2);

    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;

    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return e.id;

    size_t new_id = cell_id++;
    block_4.push_back({{v[0], v[1], v[2], v[3]}, p_dof, new_id});
    count_4++;
    return new_id;
}










Hash_ID Hash_Node::exist(size_t p0, size_t p_dof)
{
    size_t slot = hash_1(p0, p_dof) & mask_1;
    std::vector<Vertex_1>& block_1 = table_1[slot].entry_1;
    for (Vertex_1& e : block_1)
        if (e.v[0] == p0 && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}

Hash_ID Hash_Edge::exist(size_t p0, size_t p1, size_t p_dof) 
{
    if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
    size_t slot = hash_2(p0, p1, p_dof) & mask_2;      
    std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block
    for (Vertex_2& e : block_2)                               // linear scan
        if (e.v[0] == p0 && e.v[1] == p1  && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}

Hash_ID Hash_Face::exist(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;
    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}

Hash_ID Hash_Cell::exist(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;
    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}

Hash_ID Hash_Face::exist(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;
    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}

Hash_ID Hash_Cell::exist(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;
    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return {e.id, true};
    return {SIZE_MAX, false};
}









void Hash_Node::clear_table()
{
    table_1.clear();
    count_1 = 0;
    node_id = 0;
}

void Hash_Edge::clear_table()
{
    table_2.clear();
    count_2 = 0;
    edge_id = 0;
}

void Hash_Face::clear_table()
{
    table_3.clear();
    count_3 = 0;
    table_4.clear();
    count_4 = 0;
    face_id = 0;
}

void Hash_Cell::clear_table()
{
    table_3.clear();
    count_3 = 0;
    table_4.clear();
    count_4 = 0;
    cell_id = 0;
}

