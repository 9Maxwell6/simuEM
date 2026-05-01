#include "math/fem/fem_system.h"

using namespace simu;



FEM_System::FEM_System(Mesh& mesh):mesh_(mesh), block_rack_()
{
    dim_ = mesh_.get_mesh_dimension();

    block_id_ = 0;

    dof_offset_ = 0;
    dof_space_offset_ = 0;
    
    /*
    transform_triangle_o1_2D.set_mesh(mesh);
    transform_triangle_o1_3D.set_mesh(mesh);
    transform_tetrahedron_o1_3D.set_mesh(mesh);
    transform_element_2D.set_mesh(mesh);
    transform_element_3D.set_mesh(mesh);
    */

};


/**
 * @brief create hash table for volumn dof.
 * 
 * @param block provide elements with volume dof.
 * @return hash table of the volume dof based on element vertices and local dof index.
 */
const util::Block_Hash_D FEM_System::create_volume_dof_hash(const Block& block) const
{
    const std::vector<Element*>& elements = (fe_block_key_.find(block) != fe_block_key_.end()) 
                                                ? mesh_.get_element_group(fe_block_key_.at(block))
                                                : mesh_.get_mesh_elements();

    const FEM_Space * fe_space = get_block_space(block);


    util::Block_Hash_D bh_volume(32*1024);

    for (auto* e : elements) {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space = fe_space->get_basis_space(to_basis_shape(g));

        int n_dof_per_volume = fe_basis_space->get_n_dof_per_volume();

        for (int j = 0; j < n_dof_per_volume; ++j) bh_volume.get_id(idx, n, j);        
    }

    return bh_volume;
}


/**
 * @brief generate dof table given the function space and element group.
 * if using default key {0, 0}, then global elements will be used.
 * 
 * each node/edge/face/volume will be assigned to a block dof based of the order:
 *      node -> edge -> face -> volume.
 * 
 * however this information will not be stored explicitly, instead, we initialize
 * std::vector<size_t> with the internal structure:
 *      [{dof of element 1}, {dof of element 2}, {dof of element 3}, ..., {dof of element n}]
 * 
 * @example if element 1 is tetrahedron, with 4 node, 6 edge, 4 face, and 
 * 1 dof per node, 2 dof per edge, 1 dof per face, 3 dof in interior(volume). Then 
 * {dof of element 1} will be initialized as:
 * {n_1, n_2, n_3, n_4, e_1_1, e_1_2, e_2_1, e_2_2, e_3_1, e_3_2, f_1, f_2, f_3, f_4, v_1, v_2, v_3}
 * where unique dof within the block will be assigned to every components. 
 * The index directly linked to the position in block matrix
 * 
 *         [ [nodes]    .       .        .     ]
 *         [    .    [edges]    .        .     ]
 *         [    .       .    [faces]     .     ]
 *         [    .       .       .    [volumes] ]
 * 
 * notice that block index start from {0} to the {total number of dof - 1}, rather than the global index.
 * global index will be determined when assemble mutiple block matrix into the final global matrix.
 * 
 * then store this std::vector<size_t> as value of fe_block_dof_, and with block as key.
 * 
 * 
 *
 * @param block symbol of matrix block, used as key to fetch function space and mesh elements.
 * 
 * @return true if dof for this block is succeffully initialized.
 */
bool FEM_System::generate_block_dof(Block& block)
{
    FEM_Space * fe_space;
    auto it = fe_block_space_.find(block);
    if (it != fe_block_space_.end()){
        fe_space = it->second;
    }else{
        Logger::error("FEM_System::generate_block_dof - block not found: return false");
        return false;
    }

    #ifdef LOAD_PETSC
    Logger::warning("FEM_System::generate_block_dof - hash map use size_t index, while dof_idx use PetscInt.");
    #endif
     
    //FEM_Space * fe_space = fe_block_space_.at(block);

    const std::vector<Element*>& elements = (fe_block_key_.find(block) != fe_block_key_.end()) 
                                                ? mesh_.get_element_group(fe_block_key_.at(block))
                                                : mesh_.get_mesh_elements();
    
    const std::vector<Basis_Shape>& basis_shapes = fe_space->get_basis_shapes();

    bool is_node_dof   = false;
    bool is_edge_dof   = false;
    bool is_face_dof   = false;
    bool is_volume_dof = false;

    // used for initialization of hash table size
    size_t initial_size_1 = 0;
    size_t initial_size_2 = 0;
    size_t initial_size_3 = 0;
    size_t initial_size_4 = 0;


    for(Basis_Shape s : basis_shapes)
    {
        FEM_Space * fe_basis_space = fe_space->get_basis_space(s);

        int n_dof_per_node = fe_basis_space->get_n_dof_per_node();
        int n_dof_per_edge = fe_basis_space->get_n_dof_per_edge();
        int n_dof_per_face = fe_basis_space->get_n_dof_per_face();
        int n_dof_per_volume = fe_basis_space->get_n_dof_per_volume();

        if (n_dof_per_node > 0)   is_node_dof = true;
        if (n_dof_per_edge > 0)   is_edge_dof = true;
        if (n_dof_per_face > 0)   is_face_dof = true;
        if (n_dof_per_volume > 0) is_volume_dof = true;
    }

    if(is_node_dof) initial_size_1 = 32*1024;  // initialize hash table size for node
    if(is_edge_dof) initial_size_2 = 32*1024;  // initialize hash table size for edge
    if(is_face_dof) initial_size_3 = 32*1024;  // initialize hash table size for face


    // block hash table is only used for node/edge/face
    // no need for volume since there is no intersection between elements.
    util::Block_Hash bh(initial_size_1, initial_size_2, initial_size_3, initial_size_4);


    size_t element_dof_list_counter = 0; // total size of element dof list to be construct.

    bool shape_found_flag = true;

    // dof handler for creating hash table
    auto dof_hash_table_handler = [&](size_t offset, auto... args) { bh.get_id(args...); };
    

    // construct dof hash table for node/edge/face/volume separately
    // if single node/edge/face/volume has multiple dof,
    // we will store same dof id for those dof.
    // in the final step of constructing block dof list, 
    // offset will be applied for edge/face/volume dof index.
    for (auto* e : elements) 
    {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space = fe_space->get_basis_space(to_basis_shape(g));

        int n_dof_per_node   = fe_basis_space->get_n_dof_per_node();
        int n_dof_per_edge   = fe_basis_space->get_n_dof_per_edge();
        int n_dof_per_face   = fe_basis_space->get_n_dof_per_face();
        int n_dof_per_volume = fe_basis_space->get_n_dof_per_volume();

        element_dof_list_counter += fe_basis_space->get_n_node()*n_dof_per_node + 
                                    fe_basis_space->get_n_edge()*n_dof_per_edge + 
                                    fe_basis_space->get_n_face()*n_dof_per_face + 
                                    fe_basis_space->get_n_volume()*n_dof_per_volume;

        shape_found_flag = node_edge_face_dof_handler(dof_hash_table_handler, to_basis_shape(g), idx, n, 
                                                                n_dof_per_node,  0,
                                                                n_dof_per_edge,  0, 
                                                                n_dof_per_face,  0);          
    }

    // construction of hash table is now complete.


    

    // get offset for edge/face/volume dof, node dof index start from 0 hence does not need offset.
    size_t offset_edge   = bh.get_node_count();
    size_t offset_face   = bh.get_node_count() + bh.get_edge_count();
    size_t offset_volume = bh.get_node_count() + bh.get_edge_count() + bh.get_face_count();

    std::vector<dof_idx> fe_block_dof;
    fe_block_dof.reserve(element_dof_list_counter);

    size_t volume_dof_counter = 0; // volume always unique



    // dof handler for creating block dof using hash table
    auto get_id_handler = [&](size_t offset, auto... args) { fe_block_dof.push_back(offset + bh.get_id(args...)); };

    for (auto* e : elements) 
    {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space = fe_space->get_basis_space(to_basis_shape(g));

        int n_dof_per_node = fe_basis_space->get_n_dof_per_node();
        int n_dof_per_edge = fe_basis_space->get_n_dof_per_edge();
        int n_dof_per_face = fe_basis_space->get_n_dof_per_face();
        int n_dof_per_volume = fe_basis_space->get_n_dof_per_volume();

        shape_found_flag = node_edge_face_dof_handler(get_id_handler, to_basis_shape(g), idx, n, 
                                                n_dof_per_node,  0,
                                                n_dof_per_edge,  offset_edge, 
                                                n_dof_per_face,  offset_face);

        for (int j = 0; j < n_dof_per_volume; ++j) fe_block_dof.push_back(offset_volume + volume_dof_counter++);


    }
    bh.set_volume_count(volume_dof_counter);


    size_t dof_size = bh.get_node_count() + bh.get_edge_count() + bh.get_face_count() + bh.get_volume_count();
    block.row_size = dof_size;
    block.col_size = dof_size;
    fe_block_dof_data_[block] = std::move(fe_block_dof);
    fe_block_hash_[block] = std::move(bh);

    fe_block_dof_[block] = &fe_block_dof_data_[block];


    if(!shape_found_flag)
    {
        Logger::error("FEM_System::generate_block_dof - element group contains non-supported geometry.");
        return false;
    }


    return true;    
}


/**
 * @brief generate two dof table for the coupling block, using the hash table from each base block.
 * 
 * block_1 and block_2 should have shared elements where coupling happens, then:
 *      row dof of [block] should match with [ block_1 ]
 *      column dof of [block] should match with [ block_1 ]
 * 
 *      [block_1]    [block]
 *                  [block_2]
 * 
 *
 * @param block coupling block to be initialized.
 * @param block_1 initialized base block 1.
 * @param block_2 initialized base block 2.
 * 
 * @return true if dof for this block is succeffully initialized.
 */
bool FEM_System::generate_coupling_block_dof(Block& block)
{
    const std::array<const Block*, 2>& block_list = get_coupled_block(block);
    const Block& block_1 = *block_list[0];
    const Block& block_2 = *block_list[1];

    if((!block_1.is_base_block) || (!block_2.is_base_block))
    {
        Logger::error("FEM_System::generate_coupling_block_dof - block_1 or block_2 not base block, hence no dof hash table.");
        return false;
    }
    
    const std::vector<Element*>& elements = (fe_block_key_.find(block) != fe_block_key_.end()) 
                                                ? mesh_.get_element_group(fe_block_key_.at(block))
                                                : mesh_.get_mesh_elements();

    
    const util::Block_Hash& bh_1 = get_block_hash(block_1);
    const util::Block_Hash& bh_2 = get_block_hash(block_2);

    size_t n_node_dof_1   = bh_1.get_node_count();
    size_t n_edge_dof_1   = bh_1.get_edge_count();
    size_t n_face_dof_1   = bh_1.get_face_count();
    size_t n_volume_dof_1 = bh_1.get_volume_count();

    size_t n_node_dof_2   = bh_2.get_node_count();
    size_t n_edge_dof_2   = bh_2.get_edge_count();
    size_t n_face_dof_2   = bh_2.get_face_count();
    size_t n_volume_dof_2 = bh_2.get_volume_count();

    // get offset for edge/face/volume dof, node dof index start from 0 hence does not need offset.
    size_t offset_edge_1   = bh_1.get_node_count();
    size_t offset_face_1   = bh_1.get_node_count() + bh_1.get_edge_count();
    size_t offset_volume_1 = bh_1.get_node_count() + bh_1.get_edge_count() + bh_1.get_face_count();

    size_t offset_edge_2   = bh_2.get_node_count();
    size_t offset_face_2   = bh_2.get_node_count() + bh_2.get_edge_count();
    size_t offset_volume_2 = bh_2.get_node_count() + bh_2.get_edge_count() + bh_2.get_face_count();

    
    // hash table for volume
    util::Block_Hash_D bh_volume_1;
    util::Block_Hash_D bh_volume_2;
    if(n_volume_dof_1 > 0) bh_volume_1 = create_volume_dof_hash(block_1);
    if(n_volume_dof_2 > 0) bh_volume_2 = create_volume_dof_hash(block_2);
    




    const std::array<const FEM_Space *, 2>& fe_space_list = get_coupled_block_space(block);
    const FEM_Space * fe_space_1 = fe_space_list[0];
    const FEM_Space * fe_space_2 = fe_space_list[1];

    const std::vector<dof_idx>* fe_block_dof_1 = get_block_dof(block_1);
    const std::vector<dof_idx>* fe_block_dof_2 = get_block_dof(block_2);


    std::vector<dof_idx> fe_shared_block_dof_1;
    std::vector<dof_idx> fe_shared_block_dof_2;

    fe_shared_block_dof_1.reserve(fe_block_dof_1->size());
    fe_shared_block_dof_2.reserve(fe_block_dof_2->size());

    bool error_dof_flag = false;
    bool shape_found_flag = true;

    // lambda function for cheking if element ids are actually contained in hash 
    // for block hash 1
    auto get_id_handler_1 = [&](size_t offset, auto... args) {
        if (std::optional<size_t> id_o = bh_1.get_exist_id(args...)) {
            fe_shared_block_dof_1.push_back(offset + *id_o);
        } else {
            error_dof_flag = true;
        }
    };

    // for block hash 1
    auto get_id_handler_2 = [&](size_t offset, auto... args) {
        if (std::optional<size_t> id_o = bh_2.get_exist_id(args...)) {
            fe_shared_block_dof_2.push_back(offset + *id_o);
        } else {
            error_dof_flag = true;
        }
    };

    
    for (auto* e : elements) {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space_1 = fe_space_1->get_basis_space(to_basis_shape(g));
        FEM_Space * fe_basis_space_2 = fe_space_2->get_basis_space(to_basis_shape(g));

        int n_dof_per_node_1 = fe_basis_space_1->get_n_dof_per_node();
        int n_dof_per_edge_1 = fe_basis_space_1->get_n_dof_per_edge();
        int n_dof_per_face_1 = fe_basis_space_1->get_n_dof_per_face();
        int n_dof_per_volume_1 = fe_basis_space_1->get_n_dof_per_volume();

        int n_dof_per_node_2 = fe_basis_space_2->get_n_dof_per_node();
        int n_dof_per_edge_2 = fe_basis_space_2->get_n_dof_per_edge();
        int n_dof_per_face_2 = fe_basis_space_2->get_n_dof_per_face();
        int n_dof_per_volume_2 = fe_basis_space_2->get_n_dof_per_volume();

        shape_found_flag = node_edge_face_dof_handler(get_id_handler_1, to_basis_shape(g), idx, n, 
                                            n_dof_per_node_1,  0,
                                            n_dof_per_edge_1,  offset_edge_1, 
                                            n_dof_per_face_1,  offset_face_1);

        shape_found_flag = node_edge_face_dof_handler(get_id_handler_2, to_basis_shape(g), idx, n, 
                                            n_dof_per_node_2,  0,
                                            n_dof_per_edge_2,  offset_edge_2, 
                                            n_dof_per_face_2,  offset_face_2);

        for (int j = 0; j < n_dof_per_volume_1; ++j)
        {
            if (std::optional<size_t> id_o = bh_volume_1.get_exist_id(idx, n, j)){
                fe_shared_block_dof_1.push_back(offset_volume_1 + *id_o);
            }else{
                error_dof_flag = true;
            }
        } 

        for (int j = 0; j < n_dof_per_volume_2; ++j)
        {
            if (std::optional<size_t> id_o = bh_volume_2.get_exist_id(idx, n, j)){
                fe_shared_block_dof_2.push_back(offset_volume_2 + *id_o);
            }else{
                error_dof_flag = true;
            }
        }

      
    }

    fe_shared_block_dof_1.shrink_to_fit();
    fe_shared_block_dof_2.shrink_to_fit();

    if(error_dof_flag)
    {
        Logger::error("FEM_System::generate_coupling_block_dof - element group provided by the key contain elements not shared by both base blocks.");
        return false;
    }

    if(!shape_found_flag)
    {
        Logger::error("FEM_System::generate_coupling_block_dof - element group contains non-supported geometry.");
        return false;
    }


    block.row_size = block_1.row_size;
    block.col_size = block_2.col_size;

    // store actual dof list
    std::array<std::vector<dof_idx>,2>& coupled_block_dof_data = coupled_block_dof_data_[block];
    coupled_block_dof_data[0] = std::move(fe_shared_block_dof_1);
    coupled_block_dof_data[1] = std::move(fe_shared_block_dof_2);

    std::array<const std::vector<dof_idx> * ,2>& coupled_block_dof = coupled_block_dof_[block];
    coupled_block_dof[0] = &coupled_block_dof_data[0];
    coupled_block_dof[1] = &coupled_block_dof_data[1];

    /*
    size_t i=0;
    size_t j=0;
    for (auto* e : elements) {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();
        std::cout<< idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << std::endl;
        std::cout<< (*coupled_block_dof[0])[j] << " " << (*coupled_block_dof[0])[j+1] << " " << (*coupled_block_dof[0])[j+2] << " " << (*coupled_block_dof[0])[j+3] << " " << (*coupled_block_dof[0])[j+4] << " " << (*coupled_block_dof[0])[j+5]  <<std::endl;
        std::cout<< (*coupled_block_dof[1])[i] << " " << (*coupled_block_dof[1])[i+1] << " " << (*coupled_block_dof[1])[i+2] << " " << (*coupled_block_dof[1])[i+3]  <<std::endl;
        std::cout<< "-----------------------------------" << std::endl;
        i+=4;
        j+=6;
    }
    std::cout<< elements.size() << std::endl;
    //*/

    return true;

}



/**
 * @brief assign functional space to specific group of elements, 
 * if using default key {0, 0}, then assign space to global domain.
 *
 * @param fe_space finite element space.
 * @param group_key group key.
 * 
 * @return partially initialized block (with unique id, row_size and col_size, 
 * but row_offset and col_offset set to zero).
 */
Block FEM_System::register_FE_space(FEM_Space& fe_space, const Key group_key)
{
    for(const auto& [type, size] : mesh_.get_mesh_element_geometry_size())
    {
        fe_space.add_basis_shape(to_basis_shape(type));
    }

    // create uninitialized block
    block_id_++;
    Block new_block = {block_id_, 0, 0, 0, 0, true};

    fe_block_space_[new_block] = &fe_space;
    
    if(group_key.dim == 0 && group_key.id==0){
        global_space_.push_back(&fe_space);
    }else{
        group_space_[group_key].push_back(&fe_space);
        fe_block_key_[new_block] = group_key;
    }


    if(generate_block_dof(new_block)){
        return new_block;
    }else{
        Logger::error("FEM_System::register_FE_space - block initialization failed, return bad block.");
        return {0,0,0,0,0,true};
    }
}

/**
 * @brief given two group of elements under separate function space,
 * determine the shared dof index from each group. The resulting block matrix
 * will be the coupling between two function space.
 * 
 * @param block_1 block 1 with function space and group of elements.
 * @param block_2 block 1 with function space and group of elements.
 * @param shared_group_key mesh_key to group of elements that is shared by both blocks,
 * 
 * @note it is your duty to make sure the group of elements provided by shared_group_key
 * is actually shared by both blocks! If key not provided, this function will assume the
 * coupling happens across the whole mesh domain.
 * 
 * @return partially initialized block (with unique id, row_size and col_size, 
 * but row_offset and col_offset set to zero).
 */
Block FEM_System::register_FE_space_coupling(const Block& block_1, const Block& block_2, const Key shared_group_key)
{
    FEM_Space * fe_space_1;
    FEM_Space * fe_space_2;
    auto it_1 = fe_block_space_.find(block_1);
    auto it_2 = fe_block_space_.find(block_2);
    if (it_1 != fe_block_space_.end() && it_2 != fe_block_space_.end()){
        fe_space_1 = it_1->second;
        fe_space_2 = it_2->second;
    }else{
        Logger::error("FEM_System::register_FE_space_coupling - block not found: return bad block");
        return {0,0,0,0,0,false};
    }

    // create uninitialized block
    block_id_++;
    Block new_block = {block_id_, 0, 0, 0, 0, false};

    std::array<const Block*, 2>& block_pair =  coupled_block_[new_block];
    block_pair[0] = &block_1;
    block_pair[1] = &block_2;

    std::array<const FEM_Space *, 2>& block_space_pair = coupled_block_space_[new_block];
    block_space_pair[0] = fe_space_1;
    block_space_pair[1] = fe_space_2;

    if(shared_group_key.dim != 0 && shared_group_key.id!=0){
        fe_block_key_[new_block] = shared_group_key;
    }

    if(generate_coupling_block_dof(new_block)){
        return new_block;
    }else{
        Logger::error("FEM_System::register_FE_space_coupling - coupling block initialization failed, return bad block.");
        return {0,0,0,0,0,false};
    }    
}


/**
 * @brief apply Dirichlet BC to the target block.
 * 
 * @param block target block where Dirichlet BC will apply to.
 * @param group_key mesh_key to group of elements, which are the boundary of the element group in target block.
 * @param bc_type type of Dirichlet BC - Homogeneous, Constant, Field. 
 * 
 * @note this function require the hash table for this block, do not use FEM_System::delete_block_hash()
 * before using this function.
 * 
 * @return Dirichlet_BC: the struct contain all info needed to apply Dirichlet BC on the given block matrix.
 */
Dirichlet_BC FEM_System::register_Dirichlet_BC(const Block& block, const Key& group_key, Dirichlet_Type bc_type)
{
    const std::vector<Element*>& bdr_elements      = mesh_.get_element_group(group_key);
    const std::map<Geometry, size_t>& geometry_map = mesh_.get_element_geometry_size_group(group_key);
    
    

    FEM_Space* fe_space = const_cast<FEM_Space*>(get_block_space(block));
    for (const auto& [geometry, count] : geometry_map) fe_space->add_basis_shape(to_basis_shape(geometry));

    const util::Block_Hash& bh = get_block_hash(block);

    // get offset for edge/face/volume dof, node dof index start from 0 hence does not need offset.
    size_t offset_edge   = bh.get_node_count();
    size_t offset_face   = bh.get_node_count() + bh.get_edge_count();

    std::vector<dof_idx> bc_element_global_dof;

    bool error_dof_flag = false;

    // lambda function for cheking if element ids are actually contained in hash 
    // for block hash 1
    auto get_id_handler = [&](size_t offset, auto... args) {
        if (std::optional<size_t> id_o = bh.get_exist_id(args...)) {
            bc_element_global_dof.push_back(offset + *id_o);
        } else {
            error_dof_flag = true;
        }
    };


    for (auto* e : bdr_elements) 
    {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space = fe_space->get_basis_space(to_basis_shape(g));

        int n_dof_per_node = fe_basis_space->get_n_dof_per_node();
        int n_dof_per_edge = fe_basis_space->get_n_dof_per_edge();
        int n_dof_per_face = fe_basis_space->get_n_dof_per_face();

        node_edge_face_dof_handler(get_id_handler, to_basis_shape(g), idx, n, 
                                    n_dof_per_node,  0,
                                    n_dof_per_edge,  offset_edge, 
                                    n_dof_per_face,  offset_face);
    }

    Dirichlet_BC bc{.bc_type=bc_type, .block=&block, .fe_space=fe_space, .bdr_elements=&bdr_elements};
    

    std::vector<dof_idx>& bc_dofs = bc.bc_dofs;
    std::vector<dof_idx>& bc_element_dof = bc.bc_element_dof;
    // map between global block dof index and Dirichlet BC index.
    std::unordered_map<dof_idx, dof_idx> global_to_local_index;
    global_to_local_index.reserve(bc_element_global_dof.size());

    for (dof_idx i : bc_element_global_dof) 
    {
        auto [it, inserted] = global_to_local_index.try_emplace(i, bc_dofs.size());
        if(inserted) bc_dofs.push_back(i);
        if(bc_type == Dirichlet_Type::FIELD) bc_element_dof.push_back(it->second);
    }

    if(error_dof_flag)
        Logger::error("FEM_System::register_Dirichlet_BC - boundary element group provided by the key contain elements not shared by the provided block.");
    
    
    return bc;


}


const FEM_Space* FEM_System::get_block_space(const Block& block) const
{
    auto it = fe_block_space_.find(block);
    if (it != fe_block_space_.end()) return it->second;

    Logger::error("FEM_System::get_block_space - failed: block not found, return null pointer.");
    return nullptr;
}

const Key FEM_System::get_block_group_key(const Block& block) const
{
    auto it = fe_block_key_.find(block);
    if (it != fe_block_key_.end()) return it->second;

    Logger::error("FEM_System::get_block_group_key - failed: block not found, return bad key.");
    static const Key empty = {0,0};
    return empty;
}

const std::vector<dof_idx>* FEM_System::get_block_dof(const Block& block) const
{
    auto it = fe_block_dof_.find(block);
    if (it != fe_block_dof_.end()) return it->second;

    Logger::error("FEM_System::get_block_dof - failed: block not found, return nullptr.");
    return nullptr;
}

const util::Block_Hash& FEM_System::get_block_hash(const Block& block) const
{
    auto it = fe_block_hash_.find(block);
    if (it != fe_block_hash_.end()) return it->second;

    Logger::error("FEM_System::get_block_hash - failed: block not found, return empty block hash.");
    static const util::Block_Hash empty;
    return empty;
}

const std::array<const Block*, 2>& FEM_System::get_coupled_block(const Block& block) const
{
    auto it = coupled_block_.find(block);
    if (it != coupled_block_.end()) return it->second;

    Logger::error("FEM_System::get_coupled_block - failed: block not found, return empty list of block.");
    static const std::array<const Block*, 2> empty{};
    return empty;
}

const std::array<const FEM_Space *, 2>& FEM_System::get_coupled_block_space(const Block& block) const
{
    auto it = coupled_block_space_.find(block);
    if (it != coupled_block_space_.end()) return it->second;

    Logger::error("FEM_System::get_coupled_block_space - failed: block not found, return empty list of fe_space.");
    static const std::array<const FEM_Space *, 2> empty{};
    return empty;
}


const std::array<const std::vector<dof_idx> * ,2>& FEM_System::get_coupled_block_dof(const Block& block) const
{
    auto it = coupled_block_dof_.find(block);
    if (it != coupled_block_dof_.end()) return it->second;

    Logger::error("FEM_System::get_coupled_block_dof - failed: block not found, return empty list of dof.");
    static const std::array<const std::vector<dof_idx> * ,2> empty{};
    return empty;
}



const FEM_Space* FEM_System::get_block_row_space(const Block& block) const
{
    if(block.is_base_block){ return get_block_space(block); }
    else                   { return get_coupled_block_space(block)[0]; }
}

const FEM_Space* FEM_System::get_block_col_space(const Block& block) const
{
    if(block.is_base_block){ return get_block_space(block); }
    else                   { return get_coupled_block_space(block)[1]; }
}

const std::vector<dof_idx>* FEM_System::get_block_row_dof(const Block& block) const
{
    if(block.is_base_block){ return get_block_dof(block); }
    else                   { return get_coupled_block_dof(block)[0]; }
}

const std::vector<dof_idx>* FEM_System::get_block_col_dof(const Block& block) const
{
    if(block.is_base_block){ return get_block_dof(block); }
    else                   { return get_coupled_block_dof(block)[1]; }
}


Block FEM_System::transpose_block(const Block& block)
{

    block_id_++;
    Block new_block = {block_id_, 0, 0, block.col_size, block.row_size, block.is_base_block};


    if(!block.is_base_block)
    {
        const std::array<const Block*, 2>& block_list = get_coupled_block(block);
        const Block& block_1 = *block_list[0];
        const Block& block_2 = *block_list[1];

        std::array<const Block*, 2>& block_pair =  coupled_block_[new_block];
        block_pair[0] = &block_2;
        block_pair[1] = &block_1;

        const std::array<const FEM_Space *, 2>& fe_space_list = get_coupled_block_space(block);
        const FEM_Space * fe_space_1 = fe_space_list[0];
        const FEM_Space * fe_space_2 = fe_space_list[1];

        std::array<const FEM_Space *, 2>& block_space_pair = coupled_block_space_[new_block];
        block_space_pair[0] = fe_space_2;
        block_space_pair[1] = fe_space_1;

        const std::array<const std::vector<dof_idx> * ,2>& block_dof_list = get_coupled_block_dof(block);
        const std::vector<dof_idx>* block_dof_1 = block_dof_list[0];
        const std::vector<dof_idx>* block_dof_2 = block_dof_list[1];

        std::array<const std::vector<dof_idx> *, 2>& block_dof_pair = coupled_block_dof_[new_block];
        block_dof_pair[0] = block_dof_2;
        block_dof_pair[1] = block_dof_1;

    }else{
        // unusual case
        fe_block_space_[new_block] = const_cast<FEM_Space*>(get_block_space(block));
        fe_block_dof_[new_block] = get_block_dof(block);
    }

    fe_block_key_[new_block] = get_block_group_key(block);

    return new_block;
}


/**
 * @brief dof handler for node, edge and face, it will apply dof_handler operation on every
 * node, edge and face operation.
 * 
 * example:
 * 
 *      1. lambda function: dof_hash_table_handler 
 *              - initialize dof hash table for each node, edge and face.
 * 
 *      2. lambda function: get_id_handler 
 *              - get dof index for each node, edge and face.
 * 
 * @note we don't need dof handler for volume since volume is not shared between two elements.
 * 
 * @param dof_handler lambda function of operations on each dof index.
 * @param shape shape of the basis (e.g., Tetrahedron).
 * @param node_idx list of node index in mesh.
 * @param node_size size of node_idx.
 * @param n_dof_per_node number of dof per node.
 * @param n_dof_per_edge number of dof per edge.
 * @param n_dof_per_face number of dof per face.
 * @param offset_node offset for all node dof. 
 * @param offset_edge offset for all edge dof. 
 * @param offset_face offset for all face dof. 
 * 
 * @return true if handler is successfully operated.
 */
template <typename Get_dof>
bool FEM_System::node_edge_face_dof_handler(Get_dof&& dof_handler, Basis_Shape shape, const size_t* node_idx, int node_size, 
                                                                                    int n_dof_per_node, size_t offset_node,
                                                                                    int n_dof_per_edge, size_t offset_edge, 
                                                                                    int n_dof_per_face, size_t offset_face)
{
    for (int i = 0; i < node_size; ++i) 
        for (int j = 0; j < n_dof_per_node; ++j) 
            dof_handler(offset_node, node_idx[i], j);

    switch (shape) 
    {
        case Basis_Shape::TRIANGLE:
            // initialize edge dof
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[0], node_idx[1], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[0], node_idx[2], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[1], node_idx[2], j);

            // initialize face dof
            for (int j = 0; j < n_dof_per_face; ++j) dof_handler(offset_face, node_idx[0], node_idx[1], node_idx[2], j);

            break;

        case Basis_Shape::TETRAHEDRON:
            
            // initialize edge dof
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[0], node_idx[1], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[0], node_idx[2], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[0], node_idx[3], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[1], node_idx[2], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[1], node_idx[3], j);
            for (int j = 0; j < n_dof_per_edge; ++j) dof_handler(offset_edge, node_idx[2], node_idx[3], j);

            // initialize face dof
            for (int j = 0; j < n_dof_per_face; ++j) dof_handler(offset_face, node_idx[0], node_idx[1], node_idx[2], j);
            for (int j = 0; j < n_dof_per_face; ++j) dof_handler(offset_face, node_idx[0], node_idx[1], node_idx[3], j);
            for (int j = 0; j < n_dof_per_face; ++j) dof_handler(offset_face, node_idx[0], node_idx[2], node_idx[3], j);
            for (int j = 0; j < n_dof_per_face; ++j) dof_handler(offset_face, node_idx[1], node_idx[2], node_idx[3], j);

            break;

        default: 
            Logger::error("FEM_System::assign_dof - unknown Basis_Shape: return false");
            return false;
    }

    return true;
}


Block_Rack FEM_System::initialize_block_rack(size_t n_row, size_t n_col)
{
    return Block_Rack(n_row, n_col);
}



/**
 * Computes the number of nonzeros per row for a sparse matrix block arising
 * from a bilinear form between two FEM spaces (e.g. H1 × Hcurl). This is for
 * sparse matrix preallocation to avoid copying matrix data.
 *
 * Algorithm:
 *   1. Build an inverse map (row_dof_to_elem) from each row_dof to the elements that
 *      contain it, along with the starting offset into block_col_dof for that element.
 *   2. For each row_dof, iterate over its elements and collect unique col_dof
 *      using a marker array.
 *
 * @param space_1         FEM space for the row (test) functions.
 * @param block_row_dof   flat global row_dof indices for all elements, concatenated.
 * @param block_row_size  total number of unique row_dof in this block.
 * 
 * @param space_2         FEM space for the column (trial) functions.
 * @param block_col_dof   flat global col_dof indices for all elements, concatenated.
 * @param block_col_size  total number of unique col_dof in this block.
 * 
 * @param elements        list of mesh elements contributing to this block.
 *
 * @return vector of length block_row_size, where entry i is the number of
 *         nonzero columns in row i of the sparse matrix.
 */
std::vector<size_d> FEM_System::compute_nnz_per_row(
    const FEM_Space* space_1, const std::vector<dof_idx>* block_row_dof, size_d block_row_size, 
    const FEM_Space* space_2, const std::vector<dof_idx>* block_col_dof, size_d block_col_size,
    const std::vector<Element*>* elements) const
{

    // convert to 1-based indexing, 0 as default value.
    std::vector<std::vector<std::pair<const Element*,dof_idx>>> row_dof_to_elem(block_row_size+1);  
    dof_idx col_dof_offset = 0;
    dof_idx row_dof_offset = 0;
    for(const Element* e : *elements)
    {
        Basis_Shape b_shape = to_basis_shape(e->get_geometry());

        const FEM_Space* shape_space_1 = space_1->get_basis_space(b_shape);
        const FEM_Space* shape_space_2 = space_2->get_basis_space(b_shape);

        int e_row_size = shape_space_1->get_n_dof();
        int e_col_size = shape_space_2->get_n_dof();

        for (dof_idx i=row_dof_offset; i<row_dof_offset+e_row_size; ++i) 
            row_dof_to_elem[(*block_row_dof)[i]+1].push_back({e, col_dof_offset});


        row_dof_offset += e_row_size;
        col_dof_offset += e_col_size;
    }

    std::vector<dof_idx> marker(block_col_size, 0);
    std::vector<size_d> nnz(block_row_size);

    for(dof_idx i=1; i<=block_row_size; ++i)
    {
        size_d count = 0;
        for (auto&[e, col_start] : row_dof_to_elem[i])
        {
            Basis_Shape b_shape = to_basis_shape(e->get_geometry());
            const FEM_Space* shape_space_2 = space_2->get_basis_space(b_shape);
            size_t e_col_size = shape_space_2->get_n_dof();
            for(dof_idx j=col_start; j<col_start+e_col_size; ++j)
            {
                if(marker[(*block_col_dof)[j]] != i)
                {
                    marker[(*block_col_dof)[j]] = i;
                    count++;
                }
            }
        }
        // convert back to 0-based indexing
        nnz[i-1] = count;
    }

    return nnz;
}


/**
 * @brief Create Assemble_Data structure, used for global matrix assemble.
 * 
 * @param block block matrix to be assembled.
 */
Assemble_Data FEM_System::assemble_mat_data(Block& block) 
{
    const FEM_Space* space_1;
    const FEM_Space* space_2;

    const std::vector<dof_idx>* block_row_dof;
    const std::vector<dof_idx>* block_col_dof;

    if(block.is_base_block){
        const FEM_Space* space = get_block_space(block);
        space_1 = space;
        space_2 = space;

        const std::vector<dof_idx>* block_dof = get_block_dof(block);
        block_row_dof = block_dof;
        block_col_dof = block_dof;

    }else{
        const std::array<const FEM_Space*, 2>& space_list =  get_coupled_block_space(block);
        space_1 = space_list[0];
        space_2 = space_list[1];

        const std::array<const std::vector<dof_idx>*, 2>& block_dof_list =  get_coupled_block_dof(block);
        block_row_dof = block_dof_list[0];
        block_col_dof = block_dof_list[1];
    }

    const Key& group_key = get_block_group_key(block);

    const std::vector<Element*>* elements = &mesh_.get_element_group(group_key);

    // compute number of non-zero entry per row, used for preallocate the sparse block matrix.
    std::vector<size_d> nnz = compute_nnz_per_row(space_1, block_row_dof, block.row_size, 
                                                  space_2, block_col_dof, block.col_size,
                                                  elements);
    // pre-allocate block matrix
    la_kernel::init_mat(block.row_size, block.col_size, nnz, block.mat);
    //#ifdef LOAD_PETSC
    //    petsc_util::init_petsc_matrix(block.row_size, block.col_size, nnz, block.mat);
    //#else
    //    block.mat = std::make_shared<Eigen::SparseMatrix<double>>(block.row_size, block.col_size);
    //    Eigen::VectorXi nnz_eigen(nnz.begin(), nnz.end());
    //    block.mat->reserve(nnz_eigen);
    //#endif

    return Assemble_Data{
        .mesh_dim = dim_, 
        .element_dim = static_cast<int>(group_key.dim),
        .row_size = block.row_size,
        .col_size = block.col_size,
        .mesh = &mesh_,
        .space_1 = space_1,
        .space_2 = space_2,
        .elements = elements,
        .row_dof = block_row_dof,
        .col_dof = block_col_dof,
        .block_matrix = block.mat
    };

}



/**
 * @brief Create Assemble_Data structure, used for rhs vector assemble.
 * 
 * @param block block with elements to be assembled.
 */
Assemble_Data FEM_System::assemble_vec_data(Block& block) 
{
    if(!block.is_base_block){
        Logger::error("FEM_System::assemble_vec_data - block must be base block (diagonal block).");
        return Assemble_Data{};
    }
    const FEM_Space* space_1;

    const std::vector<dof_idx>* block_row_dof;

    const FEM_Space* space = get_block_space(block);
    space_1 = space;

    const std::vector<dof_idx>* block_dof = get_block_dof(block);
    block_row_dof = block_dof;

    const Key& group_key = get_block_group_key(block);

    const std::vector<Element*>* elements = &mesh_.get_element_group(group_key);

    // pre-allocate block vector
    la_kernel::init_vec(block.row_size, block.vec);
    //#ifdef LOAD_PETSC
    //    petsc_util::init_petsc_vector(block.row_size, block.vec);
    //#else
    //    block.vec = std::make_shared<Eigen::VectorXd>(Eigen::VectorXd::Zero(block.row_size));
    //#endif

    return Assemble_Data{
        .mesh_dim = dim_, 
        .element_dim = static_cast<int>(group_key.dim),
        .row_size = block.row_size,
        .mesh = &mesh_,
        .space_1 = space_1,
        .elements = elements,
        .row_dof = block_row_dof,
        .block_vector = block.vec
    };

}