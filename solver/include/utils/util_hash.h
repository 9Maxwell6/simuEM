#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

namespace util 
{


// TODO: reimplement hash with template?



// 64-bit rotate left
inline uint64_t rotl64(uint64_t x, int r) 
{
    return (x << r) | (x >> (64 - r));
}

// bit-mixer
inline size_t splitmix64(uint64_t x) 
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return static_cast<size_t>(x ^ (x >> 31));
}


struct Vertex_1 { size_t v[1]; size_t dof; size_t id; };
struct Vertex_2 { size_t v[2]; size_t dof; size_t id; };
struct Vertex_3 { size_t v[3]; size_t dof; size_t id; };
struct Vertex_4 { size_t v[4]; size_t dof; size_t id; };




struct Block_1 { std::vector<Vertex_1> entry_1; };
struct Block_2 { std::vector<Vertex_2> entry_2; };
struct Block_3 { std::vector<Vertex_3> entry_3; };
struct Block_4 { std::vector<Vertex_4> entry_4; };




// hashing constant using golden ratio
inline size_t hash_1(size_t p0, size_t p_dof=0)
{ 
    return splitmix64(p0 ^ rotl64(p_dof, 31));
}

inline size_t hash_2(size_t p0, size_t p1, size_t p_dof=0) 
{ 
    size_t h = p0 ^ rotl64(p1, 21) ^ rotl64(p_dof, 42);
    return splitmix64(h);
}

inline size_t hash_3(size_t p0, size_t p1, size_t p2, size_t p_dof=0) 
{
    size_t h = p0 ^ rotl64(p1, 16) ^ rotl64(p2, 32) ^ rotl64(p_dof, 48);
    return splitmix64(h);
}

inline size_t hash_4(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof=0) 
{ 
    size_t h = p0 ^ rotl64(p1, 10) ^ rotl64(p2, 20) ^ rotl64(p3, 30) ^ rotl64(p_dof, 50);
    return splitmix64(h);
}



struct Block_Hash 
{
    std::vector<Block_1> table_1;
    std::vector<Block_2> table_2;
    std::vector<Block_3> table_3;
    std::vector<Block_4> table_4;

   
    // mask for allocating block entry
    size_t mask_1;
    size_t mask_2;
    size_t mask_3;
    size_t mask_4;

    // increment on each insert, for generating new id and rehash
    size_t count_1 = 0;
    size_t count_2 = 0; 
    size_t count_3 = 0;
    size_t count_4 = 0;

    Block_Hash(size_t initial_size_1 = 32*1024, 
               size_t initial_size_2 = 32*1024, 
               size_t initial_size_3 = 32*1024, 
               size_t initial_size_4 = 32*1024) 
              :table_1(initial_size_1), mask_1(initial_size_1 - 1),
               table_2(initial_size_2), mask_2(initial_size_2 - 1),
               table_3(initial_size_3), mask_3(initial_size_3 - 1),
               table_4(initial_size_4), mask_4(initial_size_4 - 1)                                 
    {
        // power of 2
        assert((initial_size_1 & (initial_size_1 - 1)) == 0); 
        assert((initial_size_2 & (initial_size_2 - 1)) == 0); 
        assert((initial_size_3 & (initial_size_3 - 1)) == 0); 
        assert((initial_size_4 & (initial_size_4 - 1)) == 0); 
    }

    // 1 vertex
    size_t get_id(size_t p0, size_t p_dof)
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

    // 2 vertices
    size_t get_id(size_t p0, size_t p1, size_t p_dof) 
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

    // 3 vertices
    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p_dof)
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

    // 4 vertices
    size_t get_id(size_t p0, size_t p1, size_t p2, size_t p3, size_t p_dof)
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


    // 1 vertex
    bool if_exist(size_t p0, size_t p_dof)
    {
        size_t slot = hash_1(p0, p_dof) & mask_1;
        std::vector<Vertex_1>& block_1 = table_1[slot].entry_1;
        for (Vertex_1& e : block_1)
            if (e.v[0] == p0 && e.dof==p_dof)
                return true;
        return false;
    }

    // 2 vertices
    bool if_exist(size_t p0, size_t p1, size_t p_dof) 
    {
        if (p0 > p1) std::swap(p0, p1);                           // sort in ascending order
        size_t slot = hash_2(p0, p1, p_dof) & mask_2;      
        std::vector<Vertex_2>&  block_2 = table_2[slot].entry_2;  // get block
        for (Vertex_2& e : block_2)                               // linear scan
            if (e.v[0] == p0 && e.v[1] == p1  && e.dof==p_dof)
                return true;
        return false;
    }

    // 3 vertices
    bool if_exist(size_t p0, size_t p1, size_t p2, size_t p_dof)
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

    // 4 vertices
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



private:

    void rehash_1(size_t new_size)
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

    void rehash_2(size_t new_size)
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

    void rehash_3(size_t new_size)
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

    void rehash_4(size_t new_size)
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
    
};




template<size_t N>
struct Vertex_N { size_t v[N]; size_t dof; size_t id; };

template<size_t N>
struct Block_N { std::vector<Vertex_N<N>> entry_N; };

template <size_t N>
inline size_t hash_N(const size_t* p_list, size_t p_dof = 0) {
    uint64_t h = 0;
    for (size_t i = 0; i < N; ++i) {
        h ^= rotl64(p_list[i], (i * 17) % 64);
    }
    
    // Mix in the DOF using the next available rotation slot
    h ^= rotl64(p_dof, (N * 17) % 64);
    
    return splitmix64(h);
}


template <size_t N>
struct Block_Hash_N 
{
    std::vector<Block_N<N>> table_N;
   
    // mask for allocating block entry
    size_t mask_N;

    // increment on each insert, for generating new id and rehash
    size_t count_N = 0;

    Block_Hash_N(size_t initial_size_N = 0) :table_N(initial_size_N), mask_N(initial_size_N - 1)                                          
    {
        // power of 2
        assert((initial_size_N & (initial_size_N - 1)) == 0); 
    }

    // N vertices
    size_t get_id(size_t p[N], size_t p_dof)
    {
        if (count_N > 4 * table_N.size()) { rehash_N(table_N.size() * 2);}

        size_t v[N];
        for (size_t i = 0; i < N; ++i) v[i] = p[i];

        std::sort(v, v + N);

        size_t slot = hash_N<N>(v, p_dof) & mask_N;
        auto& block_N = table_N[slot].entry_N;

        // 5. SEARCH BUCKET
        for (auto& e : block_N) {
            // e is a Vertex_N<N>, so e.v is size_t[N]
            if (p_dof == e.dof) {
                bool match = true;
                for (size_t i = 0; i < N; ++i) {
                    if (v[i] != e.v[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return e.id;
            }
        }

        size_t new_id = count_N;
        block_N.push_back({v, p_dof, new_id});
        count_N++;
        
        return new_id;
    }

    bool if_exist(size_t p[N], size_t p_dof)
    {
        size_t v[N];
        for (size_t i = 0; i < N; ++i) {
            v[i] = p[i];
        }
        std::sort(v, v + N);
        size_t slot = hash_N<N>(v, p_dof) & mask_N;
        
        auto& block_N = table_N[slot].entry_N;

        for (const auto& e : block_N) {
            if (p_dof == e.dof) {
                bool match = true;
                for (size_t i = 0; i < N; ++i) {
                    if (e.v[i] != v[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return true;
            }
        }

        return false;
    }


private:

    void rehash_N(size_t new_size) 
    {
        std::vector<Block_N<N>> new_table(new_size);
        size_t new_mask = new_size - 1;

        for (Block_N<N>& block : table_N) {
            for (Vertex_N<N>& e : block.entry_N) {
                size_t slot = hash_N<N>(e.v, e.dof) & new_mask;
                new_table[slot].entry_N.push_back(e);
            }
        }
        table_N = std::move(new_table);
        mask_N = new_mask;
    }
    
};


// dynamic size vertex list
struct Vertex_D { std::vector<size_t> v; size_t dof; size_t id; };


struct Block_D { std::vector<Vertex_D> entry_D; };

inline size_t hash_D(const size_t* p_list, size_t N, size_t p_dof = 0) {
    uint64_t h = 0;
    for (size_t i = 0; i < N; ++i) {
        h ^= rotl64(p_list[i], (i * 17) % 64);
    }
    // Mix in the DOF based on the specific N of this element
    h ^= rotl64(p_dof, (N * 17) % 64);
    return splitmix64(h);
}

struct Block_Hash_D 
{
    std::vector<Block_D> table_D;
    size_t mask_D;
    size_t count_D = 0;

    Block_Hash_D(size_t initial_size = 32*1024) 
        : table_D(initial_size), mask_D(initial_size - 1) 
    {
        // Ensure size is a power of 2 for the bitwise mask
        assert(initial_size > 0 && (initial_size & (initial_size - 1)) == 0); 
    }

    size_t get_id(const size_t* p, size_t N, size_t p_dof)
    {
        if (count_D > 4 * table_D.size()) { rehash_D(table_D.size() * 2); }

        // 1. Local copy for sorting (preserves original mesh order)
        std::vector<size_t> v_sorted(p, p + N);
        std::sort(v_sorted.begin(), v_sorted.end());

        // 2. Locate the bucket
        size_t slot = hash_D(v_sorted.data(), N, p_dof) & mask_D;
        auto& bucket = table_D[slot].entry_D;

        // 3. Search for existing entry
        for (auto& e : bucket) {
            // Must check size N first because different element types can land in one bucket
            if (e.v.size() == N && e.dof == p_dof) {
                bool match = true;
                for (size_t i = 0; i < N; ++i) {
                    if (v_sorted[i] != e.v[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return e.id;
            }
        }

        // 4. If not found, insert new entry
        size_t new_id = count_D++;
        // Use std::move to transfer the sorted vector into the struct
        bucket.push_back({ std::move(v_sorted), p_dof, new_id });
        
        return new_id;
    }

    bool if_exist(const size_t* p, size_t N, size_t p_dof)
    {
        std::vector<size_t> v_sorted(p, p + N);
        std::sort(v_sorted.begin(), v_sorted.end());
        
        size_t slot = hash_D(v_sorted.data(), N, p_dof) & mask_D;
        auto& bucket = table_D[slot].entry_D;

        for (const auto& e : bucket) {
            if (e.v.size() == N && e.dof == p_dof) {
                bool match = true;
                for (size_t i = 0; i < N; ++i) {
                    if (e.v[i] != v_sorted[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return true;
            }
        }
        return false;
    }

private:
    void rehash_D(size_t new_size) 
    {
        std::vector<Block_D> new_table(new_size);
        size_t new_mask = new_size - 1;

        for (Block_D& old_block : table_D) {
            for (Vertex_D& e : old_block.entry_D) {
                size_t slot = hash_D(e.v.data(), e.v.size(), e.dof) & new_mask;
                new_table[slot].entry_D.push_back(std::move(e));
            }
        }
        table_D = std::move(new_table);
        mask_D = new_mask;
    }
};

}
