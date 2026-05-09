#include "math/dof/dof_manager.h"

using namespace simu;

DoF_Manager::DoF_Manager()
{
    node_h_ = Hash_Node(2);
    edge_h_ = Hash_Edge(2);
    face_h_ = Hash_Face(2);
    cell_h_ = Hash_Cell(2);
}



DoF_Manager::DoF_Manager(size_t node_size, size_t edge_size, size_t face_size, size_t cell_size)
{
    // TODO: need optimization
    node_h_ = Hash_Node(node_size);
    edge_h_ = Hash_Edge(edge_size);
    face_h_ = Hash_Face(face_size/2, face_size/2);
    cell_h_ = Hash_Cell(cell_size/2, cell_size/2);
    ready_  = true;
}

void DoF_Manager::initialize(size_t node_size, size_t edge_size, size_t face_size, size_t cell_size)
{
    // TODO: need optimization
    node_h_ = Hash_Node(node_size);
    edge_h_ = Hash_Edge(edge_size);
    face_h_ = Hash_Face(face_size/2, face_size/2);
    cell_h_ = Hash_Cell(cell_size/2, cell_size/2);
    ready_  = true;
}

void DoF_Manager::register_element_dof(Basis_Shape shape, const size_t* node_idx, 
                                       int n_dof_per_node, int n_dof_per_edge, int n_dof_per_face, int n_dof_per_cell)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        // register node dof
        for (int i = 0; i < 3; ++i) for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[i], j);
        
        // register edge dof
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j);

        // register cell dof
        for (int j = 0; j < n_dof_per_cell; ++j) cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], j);

        break;

    case Basis_Shape::TETRAHEDRON:
        // register node dof
        for (int i = 0; i < 4; ++i) for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[i], j);
        
        // register edge dof
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[3], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[3], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[2], node_idx[3], j);

        // register face dof
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], j);

        // register cell dof
        for (int j = 0; j < n_dof_per_cell; ++j) cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], node_idx[3], j);

        break;

    default: 
        Logger::error("DoF_Manager::register_element_dof - unknown Basis_Shape: return false");
    }
}










void DoF_Manager::register_node_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_node)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    for (int i = 0; i < 3; ++i) for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[i], j); break;
    case Basis_Shape::TETRAHEDRON: for (int i = 0; i < 4; ++i) for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[i], j); break;
    default: Logger::error("DoF_Manager::register_node_dof - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_edge_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_edge)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j);
        break;
    case Basis_Shape::TETRAHEDRON:
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[3], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[3], j);
        for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[2], node_idx[3], j);
        break;

    default: 
        Logger::error("DoF_Manager::register_edge_dof - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_face_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_face)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE: 
        break;

    case Basis_Shape::TETRAHEDRON:
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], j);
        for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], j);
        break;

    default: 
        Logger::error("DoF_Manager::register_face_dof - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_cell_dof(Basis_Shape shape, const size_t* node_idx, int n_dof_per_cell)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    for (int j = 0; j < n_dof_per_cell; ++j) cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], j);              break;
    case Basis_Shape::TETRAHEDRON: for (int j = 0; j < n_dof_per_cell; ++j) cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], node_idx[3], j); break;
    default: Logger::error("DoF_Manager::register_element_dof - unknown Basis_Shape: return false");
    }
}










void DoF_Manager::register_node_dof(Basis_Shape shape, const size_t* node_idx, int entity_idx, int n_dof_per_node)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    for (int i = 0; i < 3; ++i) node_h_.get_id(node_idx[i], 0); break;
        switch (entity_idx)
        {
        case 0: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[0], j); break;
        case 1: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[1], j); break;
        case 2: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[2], j); break;
        }
        break;
        
    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[0], j); break;
        case 1: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[1], j); break;
        case 2: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[2], j); break;
        case 3: for (int j = 0; j < n_dof_per_node; ++j) node_h_.get_id(node_idx[3], j); break;
        }
        break;

    default: Logger::error("DoF_Manager::register_node - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_edge_dof(Basis_Shape shape, const size_t* node_idx, int entity_idx, int n_dof_per_edge)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        switch (entity_idx)
        {
        case 0: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j); break;
        case 1: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j); break;
        case 2: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j); break;
        }
        break;
        
    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[1], j); break;
        case 1: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[2], j); break;
        case 2: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[0], node_idx[3], j); break;
        case 3: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[2], j); break;
        case 4: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[1], node_idx[3], j); break;
        case 5: for (int j = 0; j < n_dof_per_edge; ++j) edge_h_.get_id(node_idx[2], node_idx[3], j); break;
        }
        break;

    default: 
        Logger::error("DoF_Manager::register_edge - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_face_dof(Basis_Shape shape, const size_t* node_idx, int entity_idx, int n_dof_per_face)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE: 
        break;

    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], j); break;
        case 1: for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], j); break;
        case 2: for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], j); break;
        case 3: for (int j = 0; j < n_dof_per_face; ++j) face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], j); break;
        }
        break;

    default: 
        Logger::error("DoF_Manager::register_face - unknown Basis_Shape: return false");
    }
}










void DoF_Manager::register_element(Basis_Shape shape, const size_t* node_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        // register node
        for (int i = 0; i < 3; ++i) node_h_.get_id(node_idx[i], 0);
        
        // register edge
        edge_h_.get_id(node_idx[0], node_idx[1], 0);
        edge_h_.get_id(node_idx[0], node_idx[2], 0);
        edge_h_.get_id(node_idx[1], node_idx[2], 0);

        // register cell
        cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], 0);

        break;

    case Basis_Shape::TETRAHEDRON:
        // register node
        for (int i = 0; i < 4; ++i) node_h_.get_id(node_idx[i], 0);
        
        // register edge 
        edge_h_.get_id(node_idx[0], node_idx[1], 0);
        edge_h_.get_id(node_idx[0], node_idx[2], 0);
        edge_h_.get_id(node_idx[0], node_idx[3], 0);
        edge_h_.get_id(node_idx[1], node_idx[2], 0);
        edge_h_.get_id(node_idx[1], node_idx[3], 0);
        edge_h_.get_id(node_idx[2], node_idx[3], 0);

        // register face
        face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], 0);
        face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], 0);
        face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], 0);
        face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], 0);

        // register cell
        cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], node_idx[3], 0);

        break;

    default: 
        Logger::error("DoF_Manager::register_element_dof - unknown Basis_Shape: return false");
    }
}










void DoF_Manager::register_node(Basis_Shape shape, const size_t* node_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    for (int i = 0; i < 3; ++i) node_h_.get_id(node_idx[i], 0); break;
    case Basis_Shape::TETRAHEDRON: for (int i = 0; i < 4; ++i) node_h_.get_id(node_idx[i], 0); break;
    default: Logger::error("DoF_Manager::register_node - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_edge(Basis_Shape shape, const size_t* node_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        edge_h_.get_id(node_idx[0], node_idx[1], 0);
        edge_h_.get_id(node_idx[0], node_idx[2], 0);
        edge_h_.get_id(node_idx[1], node_idx[2], 0);
        break;
    case Basis_Shape::TETRAHEDRON:
        edge_h_.get_id(node_idx[0], node_idx[1], 0);
        edge_h_.get_id(node_idx[0], node_idx[2], 0);
        edge_h_.get_id(node_idx[0], node_idx[3], 0);
        edge_h_.get_id(node_idx[1], node_idx[2], 0);
        edge_h_.get_id(node_idx[1], node_idx[3], 0);
        edge_h_.get_id(node_idx[2], node_idx[3], 0);
        break;

    default: 
        Logger::error("DoF_Manager::register_edge - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_face(Basis_Shape shape, const size_t* node_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE: 
        break;

    case Basis_Shape::TETRAHEDRON:
        face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], 0);
        face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], 0);
        face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], 0);
        face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], 0);
        break;

    default: 
        Logger::error("DoF_Manager::register_face - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_cell(Basis_Shape shape, const size_t* node_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], 0);              break;
    case Basis_Shape::TETRAHEDRON: cell_h_.get_id(node_idx[0], node_idx[1], node_idx[2], node_idx[3], 0); break;
    default: Logger::error("DoF_Manager::register_cell - unknown Basis_Shape: return false");
    }
}










void DoF_Manager::register_node(Basis_Shape shape, const size_t* node_idx, int entity_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:    for (int i = 0; i < 3; ++i) node_h_.get_id(node_idx[i], 0); break;
        switch (entity_idx)
        {
        case 0: node_h_.get_id(node_idx[0], 0); break;
        case 1: node_h_.get_id(node_idx[1], 0); break;
        case 2: node_h_.get_id(node_idx[2], 0); break;
        }
        break;
        
    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: node_h_.get_id(node_idx[0], 0); break;
        case 1: node_h_.get_id(node_idx[1], 0); break;
        case 2: node_h_.get_id(node_idx[2], 0); break;
        case 3: node_h_.get_id(node_idx[3], 0); break;
        }
        break;

    default: Logger::error("DoF_Manager::register_node - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_edge(Basis_Shape shape, const size_t* node_idx, int entity_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        switch (entity_idx)
        {
        case 0: edge_h_.get_id(node_idx[0], node_idx[1], 0); break;
        case 1: edge_h_.get_id(node_idx[0], node_idx[2], 0); break;
        case 2: edge_h_.get_id(node_idx[1], node_idx[2], 0); break;
        }
        break;
        
    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: edge_h_.get_id(node_idx[0], node_idx[1], 0); break;
        case 1: edge_h_.get_id(node_idx[0], node_idx[2], 0); break;
        case 2: edge_h_.get_id(node_idx[0], node_idx[3], 0); break;
        case 3: edge_h_.get_id(node_idx[1], node_idx[2], 0); break;
        case 4: edge_h_.get_id(node_idx[1], node_idx[3], 0); break;
        case 5: edge_h_.get_id(node_idx[2], node_idx[3], 0); break;
        }
        break;

    default: 
        Logger::error("DoF_Manager::register_edge - unknown Basis_Shape: return false");
    }
}

void DoF_Manager::register_face(Basis_Shape shape, const size_t* node_idx, int entity_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE: 
        break;

    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: face_h_.get_id(node_idx[0], node_idx[1], node_idx[2], 0); break;
        case 1: face_h_.get_id(node_idx[0], node_idx[1], node_idx[3], 0); break;
        case 2: face_h_.get_id(node_idx[0], node_idx[2], node_idx[3], 0); break;
        case 3: face_h_.get_id(node_idx[1], node_idx[2], node_idx[3], 0); break;
        }
        break;

    default: 
        Logger::error("DoF_Manager::register_face - unknown Basis_Shape: return false");
    }
}










Hash_ID DoF_Manager::get_entity_id(Basis_Shape shape, const size_t* node_idx, int entity_dim, int entity_idx, int dof_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        if(entity_dim == 0){
            switch (entity_idx)
            {
            case 0: return node_h_.exist(node_idx[0], dof_idx);
            case 1: return node_h_.exist(node_idx[1], dof_idx);
            case 2: return node_h_.exist(node_idx[2], dof_idx);
            }
        }else if(entity_dim == 1){
            switch (entity_idx)
            {
            case 0: return edge_h_.exist(node_idx[0], node_idx[1], dof_idx);
            case 1: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);
            case 2: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);
            }
        }else if(entity_dim == 2){
            return cell_h_.exist(node_idx[0], node_idx[1], node_idx[2], dof_idx);
        }
        break;

    case Basis_Shape::TETRAHEDRON:
        if(entity_dim == 0){
            switch (entity_idx)
            {
            case 0: return node_h_.exist(node_idx[0], dof_idx); 
            case 1: return node_h_.exist(node_idx[1], dof_idx); 
            case 2: return node_h_.exist(node_idx[2], dof_idx); 
            case 3: return node_h_.exist(node_idx[3], dof_idx);  
            }
        }else if(entity_dim == 1){
            switch (entity_idx)
            {
            case 0: return edge_h_.exist(node_idx[0], node_idx[1], dof_idx); 
            case 1: return edge_h_.exist(node_idx[0], node_idx[2], dof_idx); 
            case 2: return edge_h_.exist(node_idx[0], node_idx[3], dof_idx);  
            case 3: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);  
            case 4: return edge_h_.exist(node_idx[1], node_idx[3], dof_idx);  
            case 5: return edge_h_.exist(node_idx[2], node_idx[3], dof_idx);  
            }
        }else if(entity_dim == 2){
            switch (entity_idx)
            {
            case 0: return face_h_.exist(node_idx[0], node_idx[1], node_idx[2], dof_idx); 
            case 1: return face_h_.exist(node_idx[0], node_idx[1], node_idx[3], dof_idx); 
            case 2: return face_h_.exist(node_idx[0], node_idx[2], node_idx[3], dof_idx);  
            case 3: return face_h_.exist(node_idx[1], node_idx[2], node_idx[3], dof_idx);  
            }
        }else if(entity_dim == 2){
            return cell_h_.exist(node_idx[0], node_idx[1], node_idx[2], node_idx[3], dof_idx);   
        }
        break;

    default: 
        Logger::error("DoF_Manager::get_entity_id - unknown Basis_Shape: return false");
    }
    return {SIZE_MAX, false};
}

Hash_ID DoF_Manager::get_node_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        switch (entity_idx)
        {
        case 0: return node_h_.exist(node_idx[0], dof_idx);
        case 1: return node_h_.exist(node_idx[1], dof_idx);
        case 2: return node_h_.exist(node_idx[2], dof_idx);
        }
        break;

    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: return node_h_.exist(node_idx[0], dof_idx); 
        case 1: return node_h_.exist(node_idx[1], dof_idx); 
        case 2: return node_h_.exist(node_idx[2], dof_idx); 
        case 3: return node_h_.exist(node_idx[3], dof_idx);  
        }
        break;

    default: 
        Logger::error("DoF_Manager::get_node_id - unknown Basis_Shape: return false");
    }
    return {SIZE_MAX, false};
}

Hash_ID DoF_Manager::get_edge_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        switch (entity_idx)
        {
        case 0: return edge_h_.exist(node_idx[0], node_idx[1], dof_idx);
        case 1: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);
        case 2: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);
        }
        break;

    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: return edge_h_.exist(node_idx[0], node_idx[1], dof_idx); 
        case 1: return edge_h_.exist(node_idx[0], node_idx[2], dof_idx); 
        case 2: return edge_h_.exist(node_idx[0], node_idx[3], dof_idx);  
        case 3: return edge_h_.exist(node_idx[1], node_idx[2], dof_idx);  
        case 4: return edge_h_.exist(node_idx[1], node_idx[3], dof_idx);  
        case 5: return edge_h_.exist(node_idx[2], node_idx[3], dof_idx);  
        }
        break;

    default: 
        Logger::error("DoF_Manager::get_edge_id - unknown Basis_Shape: return false");
    }
    return {SIZE_MAX, false};
}

Hash_ID DoF_Manager::get_face_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE:
        break;

    case Basis_Shape::TETRAHEDRON:
        switch (entity_idx)
        {
        case 0: return face_h_.exist(node_idx[0], node_idx[1], node_idx[2], dof_idx); 
        case 1: return face_h_.exist(node_idx[0], node_idx[1], node_idx[3], dof_idx); 
        case 2: return face_h_.exist(node_idx[0], node_idx[2], node_idx[3], dof_idx);  
        case 3: return face_h_.exist(node_idx[1], node_idx[2], node_idx[3], dof_idx);  
        }
        break;

    default: 
        Logger::error("DoF_Manager::get_entity_id - unknown Basis_Shape: return false");
    }
    return {SIZE_MAX, false};
}

Hash_ID DoF_Manager::get_cell_id(Basis_Shape shape, const size_t* node_idx, int entity_idx, int dof_idx)
{
    switch (shape) 
    {        
    case Basis_Shape::TRIANGLE: return cell_h_.exist(node_idx[0], node_idx[1], node_idx[2], dof_idx);
    case Basis_Shape::TETRAHEDRON: return cell_h_.exist(node_idx[0], node_idx[1], node_idx[2], node_idx[3], dof_idx);   
    default: 
        Logger::error("DoF_Manager::get_cell_id - unknown Basis_Shape: return false");
    }
    return {SIZE_MAX, false};
}










void DoF_Manager::delete_all()
{
    node_h_.clear_table();
    edge_h_.clear_table();
    face_h_.clear_table();
    cell_h_.clear_table();
}

void DoF_Manager::delete_node_hash(){ node_h_.clear_table(); }
void DoF_Manager::delete_edge_hash(){ edge_h_.clear_table(); }
void DoF_Manager::delete_face_hash(){ face_h_.clear_table(); }
void DoF_Manager::delete_cell_hash(){ cell_h_.clear_table(); }
