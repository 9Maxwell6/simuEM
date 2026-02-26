#pragma once

namespace util 
{


#include <vector>
#include <algorithm>
#include <cassert>


struct Vertex_2 { size_t v[2]; size_t id; };
struct Vertex_3 { size_t v[3]; size_t id; };
struct Vertex_4 { size_t v[4]; size_t id; };

struct Block_2 { std::vector<Vertex_2> entry_2; };
struct Block_3 { std::vector<Vertex_3> entry_3; };
struct Block_4 { std::vector<Vertex_4> entry_4; };


// hashing constant using golden ratio
inline size_t hash_2(size_t p0, size_t p1) 
    { return 0x9E3779B9ul*p0 + 0xDAA66D2Bul*p1; }

inline size_t hash_3(size_t p0, size_t p1, size_t p2) 
    { return 0x9E3779B9ul*p0 + 0xDAA66D2Bul*p1 + 0x7AD5623Dul*p2; }

inline size_t hash_4(size_t p0, size_t p1, size_t p2, size_t p3) 
    { return 0x9E3779B9ul*p0 + 0xDAA66D2Bul*p1 + 0x7AD5623Dul*p2 + 0x1904564Ful*p3; }



struct Block_Hash {
    
    std::vector<Block_2> table_2;
    std::vector<Block_3> table_3;
    std::vector<Block_4> table_4;
   
    // mask for allocating block entry
    size_t mask_2;
    size_t mask_3;
    size_t mask_4;

    // increment on each insert, for generating new id and rehash
    size_t count_2 = 0; 
    size_t count_3 = 0;
    size_t count_4 = 0;

    Block_Hash(size_t initial_size_2 = 32*1024, size_t initial_size_3 = 32*1024, size_t initial_size_4 = 32*1024) 
        : table_2(initial_size_2), mask_2(initial_size_2 - 1),
          table_3(initial_size_3), mask_3(initial_size_3 - 1),
          table_4(initial_size_4), mask_4(initial_size_4 - 1)                                          
    {
        assert((initial_size_2 & (initial_size_2 - 1)) == 0); // power of 2
        assert((initial_size_3 & (initial_size_3 - 1)) == 0); // power of 2
        assert((initial_size_4 & (initial_size_4 - 1)) == 0); // power of 2
    }

    // 2 vertices
    size_t get_id(size_t p0, size_t p1) 
    {
        // check rehash
        if (count_2 > 4*table_2.size()) rehash_2(table_2.size() * 2);

        if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
        size_t slot = hash_2(p0, p1) & mask_2;      
        std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block

        for (Vertex_2& e : block_2)                               // linear scan
            if (e.v[0] == p0 && e.v[1] == p1)
                return e.id;

        size_t new_id = count_2;                                  // global id in table
        block_2.push_back({{p0, p1}, new_id});

        count_2++;

        return new_id;
    }

    // 3 vertices
    size_t get_id(size_t p0, size_t p1, size_t p2)
    {
        if (count_3 > 4*table_3.size()) rehash_3(table_3.size() * 2);

        size_t v[3] = {p0, p1, p2};
        std::sort(v, v + 3);
        size_t slot = hash_3(v[0], v[1], v[2]) & mask_3;
        std::vector<Vertex_3>&  block_3 = table_3[slot].entry_3;

        for (Vertex_3& e : block_3)
            if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2])
                return e.id;

        size_t new_id = count_3;  
        block_3.push_back({{v[0], v[1], v[2]}, new_id});

        count_3++;

        return new_id;
    }

    // 4 vertices
    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p3)
    {
        if (count_4 > 4*table_4.size()) rehash_4(table_4.size() * 2);

        size_t v[4] = {p0, p1, p2, p3};
        std::sort(v, v + 4);
        size_t slot = hash_4(v[0], v[1], v[2], v[3]) & mask_4;
        std::vector<Vertex_4>&  block_4 = table_4[slot].entry_4;

        for (auto& e : block_4)
            if (e.v[0] == v[0] && e.v[1] == v[1] && e.v[2] == v[2] && e.v[3] == v[3])
                return e.id;

        size_t new_id = count_4;
        block_4.push_back({{v[0], v[1], v[2], v[3]}, new_id});

        count_4++;

        return new_id;
    }


private:

    void rehash_2(size_t new_size)
    {
        std::vector<Block_2> new_table(new_size);
        size_t new_mask = new_size - 1;

        for (Block_2& block : table_2)
            for (Vertex_2& e : block.entry_2) {
                size_t slot = hash_2(e.v[0], e.v[1]) & new_mask;
                new_table[slot].entry_2.push_back(e);
            }

        table_2 = std::move(new_table);
        mask_2  = new_mask;
    }

    void rehash_3(size_t new_size)
    {
        std::vector<Block_3> new_table(new_size);
        size_t new_mask = new_size - 1;

        for (Block_3& block : table_3)
            for (Vertex_3& e : block.entry_3) {
                size_t slot = hash_3(e.v[0], e.v[1], e.v[2]) & new_mask;
                new_table[slot].entry_3.push_back(e);
            }

        table_3 = std::move(new_table);
        mask_3  = new_mask;
    }

    void rehash_4(size_t new_size)
    {
        std::vector<Block_4> new_table(new_size);
        size_t new_mask = new_size - 1;

        for (Block_4& block : table_4)
            for (Vertex_4& e : block.entry_4) {
                size_t slot = hash_4(e.v[0], e.v[1], e.v[2], e.v[3]) & new_mask;
                new_table[slot].entry_4.push_back(e);
            }

        table_4 = std::move(new_table);
        mask_4  = new_mask;
    }
    
};

}