#include "utils/util_hash.h"
#include "utils/logger.h"


Block_Hash_1::Block_Hash_1(size_t initial_size = 32*1024) : table(initial_size), mask(initial_size - 1)
{
    if((initial_size & (initial_size - 1)) != 0) 
    {
        initial_size = 32*1024;
        table(initial_size);
        mask(initial_size - 1);
        Logger::warning("Block_Hash_1 - hash table size not power of 2, use default size: 32*1024");
    }
}


size_t Block_Hash_1::get_id(size_t p0, size_t p_dof)
{
    //std::cout<<p0<<std::endl;
    if (count_1 > 4*table_1.size()) rehash_1(table_1.size() * 2);

    size_t slot = hash_1(p0, p_dof) & mask_1;
    std::vector<Vertex_1>& block_1 = table_1[slot].entry_1;
    for (Vertex_1& e : block_1)
        if (e.v[0] == p0 && e.dof==p_dof)
            return e.id;

    size_t new_id = count_1;
    block_1.push_back({{p0}, p_dof, new_id});
    count_1++;
    return new_id;
}

bool Block_Hash_1::if_exist(size_t p0, size_t p_dof)
{
    size_t slot = hash_1(p0, p_dof) & mask_1;
    std::vector<Vertex_1>& block_1 = table_1[slot].entry_1;
    for (Vertex_1& e : block_1)
        if (e.v[0] == p0 && e.dof==p_dof)
            return true;
    return false;
}

void Block_Hash_1::rehash_1(size_t new_size)
{
    std::vector<Block_1> new_table(new_size);
    size_t new_mask = new_size - 1;
    for (Block_1& block : table_1)
        for (Vertex_1& e : block.entry_1) 
        {
            size_t slot = hash_1(e.v[0], e.dof) & new_mask;
            new_table[slot].entry_1.push_back(e);
        }

    table_1 = std::move(new_table);
    mask_1  = new_mask;
}

Block_Hash_2::Block_Hash_2(size_t initial_size = 32*1024) : table(initial_size), mask(initial_size - 1)
{
    if((initial_size & (initial_size - 1)) != 0) 
    {
        initial_size = 32*1024;
        table(initial_size);
        mask(initial_size - 1);
        Logger::warning("Block_Hash_2 - hash table size not power of 2, use default size: 32*1024");
    }
}

size_t Block_Hash_2::get_id(size_t p0, size_t p1, size_t p_dof) 
{
    // check rehash
    if (count_2 > 4*table_2.size()) rehash_2(table_2.size() * 2);

    if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
    size_t slot = hash_2(p0, p1, p_dof) & mask_2;      
    std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block

    for (Vertex_2& e : block_2)                               // linear scan
        if (e.v[0] == p0 && e.v[1] == p1  && e.dof==p_dof)
            return e.id;

    size_t new_id = count_2;                                  // global id in table
    block_2.push_back({{p0, p1}, p_dof, new_id});
    count_2++;
    return new_id;
}

bool Block_Hash_2::if_exist(size_t p0, size_t p1, size_t p_dof) 
{
    if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
    size_t slot = hash_2(p0, p1, p_dof) & mask_2;      
    std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block
    for (Vertex_2& e : block_2)                               // linear scan
        if (e.v[0] == p0 && e.v[1] == p1  && e.dof==p_dof)
            return true;
    return false;
}

void Block_Hash_2::rehash_2(size_t new_size)
{
    std::vector<Block_2> new_table(new_size);
    size_t new_mask = new_size - 1;

    for (Block_2& block : table_2)
        for (Vertex_2& e : block.entry_2)
        {
            size_t slot = hash_2(e.v[0], e.v[1], e.dof) & new_mask;
            new_table[slot].entry_2.push_back(e);
        }

    table_2 = std::move(new_table);
    mask_2  = new_mask;
}


Block_Hash_3::Block_Hash_3(size_t initial_size = 32*1024) : table(initial_size), mask(initial_size - 1)
{
    if((initial_size & (initial_size - 1)) != 0) 
    {
        initial_size = 32*1024;
        table(initial_size);
        mask(initial_size - 1);
        Logger::warning("Block_Hash_3 - hash table size not power of 2, use default size: 32*1024");
    }
}

size_t Block_Hash_3::get_id(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    if (count_3 > 4*table_3.size()) rehash_3(table_3.size() * 2);

    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;

    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return e.id;

    size_t new_id = count_3;  
    block_3.push_back({{v[0], v[1], v[2]}, p_dof, new_id});
    count_3++;
    return new_id;
}

bool Block_Hash_3::if_exist(size_t p0, size_t p1, size_t p2, size_t p_dof)
{
    size_t v[3] = {p0, p1, p2};
    std::sort(v, v + 3);
    size_t slot = hash_3(v[0], v[1], v[2], p_dof) & mask_3;
    std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;
    for (Vertex_3& e : block_3)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2]  && e.dof==p_dof)
            return true;
    return false;
}

void Block_Hash_3::rehash_3(size_t new_size)
{
    std::vector<Block_3> new_table(new_size);
    size_t new_mask = new_size - 1;

    for (Block_3& block : table_3)
        for (Vertex_3& e : block.entry_3)
        {
            size_t slot = hash_3(e.v[0], e.v[1], e.v[2], e.dof) & new_mask;
            new_table[slot].entry_3.push_back(e);
        }

    table_3 = std::move(new_table);
    mask_3  = new_mask;
}


Block_Hash_4::Block_Hash_4(size_t initial_size = 32*1024) : table(initial_size), mask(initial_size - 1)
{
    if((initial_size & (initial_size - 1)) != 0) 
    {
        initial_size = 32*1024;
        table(initial_size);
        mask(initial_size - 1);
        Logger::warning("Block_Hash_4 - hash table size not power of 2, use default size: 32*1024");
    }
}

size_t Block_Hash_4::get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    if (count_4 > 4*table_4.size()) rehash_4(table_4.size() * 2);

    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;

    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return e.id;

    size_t new_id = count_4;
    block_4.push_back({{v[0], v[1], v[2], v[3]}, p_dof, new_id});
    count_4++;
    return new_id;
}

void Block_Hash_4::rehash_4(size_t new_size)
{
    std::vector<Block_4> new_table(new_size);
    size_t new_mask = new_size - 1;

    for (Block_4& block : table_4)
        for (Vertex_4& e : block.entry_4)
        {
            size_t slot = hash_4(e.v[0], e.v[1], e.v[2], e.v[3], e.dof) & new_mask;
            new_table[slot].entry_4.push_back(e);
        }

    table_4 = std::move(new_table);
    mask_4  = new_mask;
}

bool if_exist(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
{
    size_t v[4] = {p0, p1, p2, p3};
    std::sort(v, v + 4);
    size_t slot = hash_4(v[0], v[1], v[2], v[3], p_dof) & mask_4;
    std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;
    for (auto& e : block_4)
        if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3]  && e.dof==p_dof)
            return true;
    return false;
}