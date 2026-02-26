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
 * @brief 
 *
 * @param fe_space finite element space.
 * @param group_key group key.
 * 
 * @return uninitialized block (with unique id, but offset/row_size/col_size set to zero).
 */
void FEM_System::generate_space_dof_table(FEM_Space& fe_space,  const Key group_key)
{
    const std::vector<Element*>& elements = (group_key.dim == 0 && group_key.id==0) ? mesh_.get_mesh_elements() : mesh_.get_element_group(group_key);
    
    const std::vector<Basis_Shape>& basis_shapes = fe_space.get_basis_shapes();

    bool is_node_dof   = false;
    bool is_edge_dof   = false;
    bool is_face_dof   = false;
    bool is_volume_dof = false;

    for(Basis_Shape s : basis_shapes)
    {
        FEM_Space * fe_basis_space = fe_space.get_basis_space(s);
        if (fe_basis_space->get_n_dof_per_node() > 0)   is_node_dof = true;
        if (fe_basis_space->get_n_dof_per_edge() > 0)   is_edge_dof = true;
        if (fe_basis_space->get_n_dof_per_face() > 0)   is_face_dof = true;
        if (fe_basis_space->get_n_dof_per_volume() > 0) is_volume_dof = true;
    }
    
    std::set<size_t>              unique_nodes;
    std::set<std::vector<size_t>> unique_edges;
    std::set<std::vector<size_t>> unique_faces;
    size_t                     num_volumes = 0; // volume always unique

    auto makeKey = [](std::initializer_list<size_t> list) {
        std::vector<size_t> v(list);
        std::sort(v.begin(), v.end());
        return v;
    };

    for (auto* e : elements) {
        const size_t* idx  = e->get_nodeIdx();
        int n = e->get_nodeNum();

        auto g = e->get_geometry();

        for (int i = 0; i < n; ++i) unique_nodes.insert(idx[i]);

        switch (g) {
            case Geometry::EDGE:
                unique_edges.insert(makeKey({idx[0], idx[1]}));
                break;

            case Geometry::TRIANGLE:
                unique_edges.insert(makeKey({idx[0], idx[1]}));
                unique_edges.insert(makeKey({idx[1], idx[2]}));
                unique_edges.insert(makeKey({idx[0], idx[2]}));
                unique_faces.insert(makeKey({idx[0], idx[1], idx[2]}));
                break;

            case Geometry::TETRAHEDRON:
                unique_edges.insert(makeKey({idx[0], idx[1]}));
                unique_edges.insert(makeKey({idx[0], idx[2]}));
                unique_edges.insert(makeKey({idx[0], idx[3]}));
                unique_edges.insert(makeKey({idx[1], idx[2]}));
                unique_edges.insert(makeKey({idx[1], idx[3]}));
                unique_edges.insert(makeKey({idx[2], idx[3]}));
                unique_faces.insert(makeKey({idx[0], idx[1], idx[2]}));
                unique_faces.insert(makeKey({idx[0], idx[1], idx[3]}));
                unique_faces.insert(makeKey({idx[0], idx[2], idx[3]}));
                unique_faces.insert(makeKey({idx[1], idx[2], idx[3]}));
                num_volumes++;
                break;


            default: break;
        }
    }

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

}