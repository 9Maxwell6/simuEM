#pragma once

#include <functional>

#include "math/dof/dof_hash.h"
#include "math/fem/shape.h"


namespace simu {


class DoF_Manager
{
private:

    Hash_Node node_h_;
    Hash_Edge edge_h_;
    Hash_Face face_h_;
    Hash_Cell cell_h_;

    bool ready_ = false;


    // sometimes you do not want to register immedietely, then you could store the operation into the buffer and execute later.
    // usage case:
    //      - for each edge, two client want to do something.
    //      - client 1 call get_edge_id, if not exist, do (...), and register_edge.
    //      - client 2 call get_edge_id, if not exist, do (...), and register_edge.
    //      - but since client 1 already registered the edge, 
    //        operation in client 2 will not be executed.
    //     -> store register_edge() into the buffer and call it later.
    std::vector<std::function<void()>> operation_buffer;


public:
    DoF_Manager();
    DoF_Manager(size_t node_size, size_t edge_size, size_t face_size, size_t cell_size);

    void initialize(size_t node_size=1024, size_t edge_size=1024, size_t face_size=1024, size_t cell_size=1024);

    bool is_ready() const { return ready_; }

    // add operation into buffer
    template<typename Func>
    void pending_operation(Func&& f)
    {
        operation_buffer.emplace_back(std::forward<Func>(f));
    }

    void execute_pending_operations()
    {
        auto ops = std::move(operation_buffer);
        operation_buffer.clear();
        for (auto& op : ops) op();
    }
   
    void register_element_dof(Basis_Shape shape, const size_t* node_idx,
                              int n_dof_per_node, int n_dof_per_edge, int n_dof_per_face, int n_dof_per_cell);
    void register_node_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_node);
    void register_edge_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_edge);
    void register_face_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_face);
    void register_cell_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_cell);

    void register_node_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_node, int entity_idx);
    void register_edge_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_edge, int entity_idx);
    void register_face_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_face, int entity_idx);

    
    // ignore dof (should not mix-use with dof-involved function)
    void register_element(Basis_Shape shape, const size_t* node_idx);
    void register_node(Basis_Shape shape, const size_t* node_idx);
    void register_edge(Basis_Shape shape, const size_t* node_idx);
    void register_face(Basis_Shape shape, const size_t* node_idx);
    void register_cell(Basis_Shape shape, const size_t* node_idx);

    void register_node(Basis_Shape shape, const size_t* node_idx, int entity_idx);
    void register_edge(Basis_Shape shape, const size_t* node_idx, int entity_idx);
    void register_face(Basis_Shape shape, const size_t* node_idx, int entity_idx);
    

    //TODO
    Hash_ID get_entity_id(Basis_Shape shape, const size_t* node_idx, int entity_dim, int entity_idx, int dof_idx=0);
    Hash_ID get_node_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx=0);
    Hash_ID get_edge_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx=0);
    Hash_ID get_face_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx=0);
    Hash_ID get_cell_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx=0);


    // delete hash table. (hash table can be very large, It is recommended to delete it as soon as it is no longer in use.)
    void delete_all();
    void delete_node_hash();
    void delete_edge_hash();
    void delete_face_hash();
    void delete_cell_hash();
};


}
