#include "math/fem/fem_system.h"

using namespace simu;



FEM_System::FEM_System(Mesh& mesh):mesh_(mesh)
{
    dim_ = mesh_.get_mesh_dimension();

    block_id_ = 0;

    dof_offset_ = 0;
    dof_space_offset_ = 0;
    

};


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
    FEM_Space * fe_space = fe_block_space_.at(block);

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

    std::vector<size_t> node_dof_id;
    std::vector<size_t> edge_dof_id;
    std::vector<size_t> face_dof_id;
    std::vector<size_t> volume_dof_id;

    // overestimate the size of final dof table (usually much larger than actually needed)
    size_t estimate_node_dof_size = 0;
    size_t estimate_edge_dof_size = 0;
    size_t estimate_face_dof_size = 0;
    size_t estimate_volume_dof_size = 0;

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
        
        estimate_node_dof_size += elements.size()*(fe_basis_space->get_n_node()*n_dof_per_node);
        estimate_edge_dof_size += elements.size()*(fe_basis_space->get_n_edge()*n_dof_per_edge);
        estimate_face_dof_size += elements.size()*(fe_basis_space->get_n_face()*n_dof_per_face);
        estimate_volume_dof_size += elements.size()*(fe_basis_space->get_n_volume()*n_dof_per_volume);
                                              

        switch (s)
        {
            case Basis_Shape::TETRAHEDRON:
                if(is_node_dof) initial_size_1 = 32*1024;  // initialize hash table size for node
                if(is_edge_dof) initial_size_2 = 32*1024;  // initialize hash table size for edge
                if(is_face_dof) initial_size_3 = 32*1024;  // initialize hash table size for face
                break;
            default:
                Logger::error("FEM_System::generate_block_dof - unknown Basis_Shape: return false");
                return false;
            
                
        }
    }

    estimate_node_dof_size /= basis_shapes.size();
    estimate_edge_dof_size /= basis_shapes.size();
    estimate_face_dof_size /= basis_shapes.size();
    estimate_volume_dof_size /= basis_shapes.size();

    std::vector<size_t> block_node_dof;
    std::vector<size_t> block_edge_dof;
    std::vector<size_t> block_face_dof;
    std::vector<size_t> block_volume_dof;

    block_node_dof.reserve(estimate_node_dof_size);
    block_edge_dof.reserve(estimate_edge_dof_size);
    block_face_dof.reserve(estimate_face_dof_size);
    block_volume_dof.reserve(estimate_volume_dof_size);


    // block hash table is only used for node/edge/face
    // no need for volume since there is no intersection between elements.
    util::Block_Hash bh(initial_size_1, initial_size_2, initial_size_3, initial_size_4);

    
    
    std::set<size_t>              unique_nodes;
    std::set<std::vector<size_t>> unique_edges;
    std::set<std::vector<size_t>> unique_faces;
    size_t                        volume_dof_counter = 0; // volume always unique

    // construct dof index for node/edge/face/volume separately
    // if single node/edge/face/volume has multiple dof,
    // we will store same dof id for those dof.
    // in the final step of constructing block dof list, 
    // those duplicate id value will be shifted until each dof has unique id.
    for (auto* e : elements) {
        const size_t* idx  = e->get_nodeIdx();
        int n = e->get_nodeNum();

        Geometry g = e->get_geometry();

        FEM_Space * fe_basis_space = fe_space->get_basis_space(to_basis_shape(g));

        int n_dof_per_node = fe_basis_space->get_n_dof_per_node();
        int n_dof_per_edge = fe_basis_space->get_n_dof_per_edge();
        int n_dof_per_face = fe_basis_space->get_n_dof_per_face();
        int n_dof_per_volume = fe_basis_space->get_n_dof_per_volume();

        // assign dof to nodes
        if(is_node_dof)
            for (int i = 0; i < n; ++i) 
                for (int j = 0; j < n_dof_per_node; ++j) 
                {
                    block_node_dof.push_back(bh.get_id(idx[i], j));
                    //std::cout<<bh.get_id(idx[i], j)<<std::endl;
                }
                //std::cout<<"------------------"<<std::endl;

        switch (to_basis_shape(g)) {
            case Basis_Shape::TETRAHEDRON:
                if(is_edge_dof)
                {
                    // TODO multiple dof per edge or face??
                    // for now just store single id, then we hand   le it in the final stage
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[0], idx[1], j));
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[0], idx[2], j));
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[0], idx[3], j));
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[1], idx[2], j));
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[1], idx[3], j));
                    for (int j = 0; j < n_dof_per_edge; ++j) block_edge_dof.push_back(bh.get_id(idx[2], idx[3], j));
                    
                }

                if(is_face_dof)
                {
                    for (int j = 0; j < n_dof_per_face; ++j) block_face_dof.push_back(bh.get_id(idx[0], idx[1], idx[2], j));
                    for (int j = 0; j < n_dof_per_face; ++j) block_face_dof.push_back(bh.get_id(idx[0], idx[1], idx[3], j));
                    for (int j = 0; j < n_dof_per_face; ++j) block_face_dof.push_back(bh.get_id(idx[0], idx[2], idx[3], j));
                    for (int j = 0; j < n_dof_per_face; ++j) block_face_dof.push_back(bh.get_id(idx[1], idx[2], idx[3], j));

                }
                
                if(is_volume_dof)
                {   
                    
                    for (int j = 0; j < n_dof_per_volume; ++j)
                    {
                        block_volume_dof.push_back(volume_dof_counter);
                        volume_dof_counter++;
                    }
                }
                break;


            default: 
                Logger::error("FEM_System::generate_block_dof - unknown Basis_Shape: return false");
                return false;
        }
    }

    block_node_dof.shrink_to_fit();
    
    block_edge_dof.shrink_to_fit();
    block_face_dof.shrink_to_fit();
    block_volume_dof.shrink_to_fit();

    std::vector<size_t> fe_block_dof;
    fe_block_dof.reserve(block_node_dof.size()+block_edge_dof.size()+block_face_dof.size()+block_volume_dof.size());
    // 2. Append each vector
    fe_block_dof.insert(fe_block_dof.end(), block_node_dof.begin(), block_node_dof.end());
    fe_block_dof.insert(fe_block_dof.end(), block_edge_dof.begin(), block_edge_dof.end());
    fe_block_dof.insert(fe_block_dof.end(), block_face_dof.begin(), block_face_dof.end());
    fe_block_dof.insert(fe_block_dof.end(), block_volume_dof.begin(), block_volume_dof.end());


    //Logger::info(std::to_string(fe_block_dof.size()));
    size_t i=0;
    for (auto* e : elements) {
        const size_t* idx  = e->get_nodeIdx();
        int n = e->get_nodeNum();
        std::cout<< idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << " | " << fe_block_dof[i] << " " << fe_block_dof[i+1] << " " << fe_block_dof[i+2] << " " << fe_block_dof[i+3]  <<std::endl;
        i+=4;
    }
    for(size_t i=0 ; i<fe_block_dof.size(); i+=4){
        //std::cout<< fe_block_dof[i] << " " << fe_block_dof[i+1] << " " << fe_block_dof[i+2] << " " << fe_block_dof[i+3]  <<std::endl;
    }
    fe_block_dof_[block] = fe_block_dof;
    return true;
    //return {unique_nodes.size(), unique_edges.size(), unique_faces.size(), num_volumes};
    
}




bool FEM_System::initialize_space_dof()
{

}



size_t FEM_System::assign_element_dof(FEM_Space& fe_space, Element& e)
{
    size_t current_offset = dof_offset_;
    FEM_Space * basis_space = fe_space.get_basis_space(to_basis_shape(e.get_geometry()));
    

    return current_offset;
}




Basis_Shape FEM_System::to_basis_shape(Geometry t)
{
    switch (t) {
        case Geometry::TETRAHEDRON: return Basis_Shape::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - type not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

Geometry FEM_System::to_element_geometry(Basis_Shape g)
{
    switch (g) {
        case Basis_Shape::TETRAHEDRON: return Geometry::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - geometry not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}


/*
bool FEM_System::create_FE_space(FEM_Space * fe_space){

    std::vector<FEM_Space *> new_fem_space;

    switch (dim_)
    {
    case 3:
        //if(mesh_.)
        break;
    case 2:
        break;
    case 1:
        break;
    default:
        Logger::error("FEM_System::create_FE_space - impossible dimension: "+std::to_string(dim_));
        
    }
}
    */


/**
 * @brief Assign space to all elements in mesh
 *
 * @param fs finite element space.
 * @param p_order order of basis function.
 */
//bool FEM_System::initialize_FE_space(FEM_Space& fe_space)
//{
    // TODO: should use vector<Space>, what if we assign duplicate Space?
    // use vector<FEM_Space *>,
    //global_space_.push_back(fs);
    
//}


/**
 * @brief assign functional space to specific group of elements, 
 * if using default key {0, 0}, then assign space to global domain.
 *
 * @param fe_space finite element space.
 * @param group_key group key.
 * 
 * @return uninitialized block (with unique id, but offset/row_size/col_size set to zero).
 */
Block FEM_System::register_FE_space(FEM_Space& fe_space, const Key group_key)
{
    for(const auto& [type, size] : mesh_.get_mesh_element_geometry_size())
    {
        fe_space.add_basis_shape(to_basis_shape(type));
    }

    // create uninitialized block
    block_id_++;
    Block new_block = {block_id_, 0, 0, 0};

    fe_block_space_[new_block] = &fe_space;
    
    if(group_key.dim == 0 && group_key.id==0){
        global_space_.push_back(&fe_space);
    }else{
        group_space_[group_key].push_back(&fe_space);
        fe_block_key_[new_block] = group_key;
    }

    generate_block_dof(new_block);

}