#pragma once

#include "fem_space.h"


namespace simu {



class Hcurl_Space : public FEM_Space 
{
private:
    //std::vector<std::unique_ptr<Hcurl_Space>> shape_Hcurl_;
    std::unordered_map<Basis_Shape, std::unique_ptr<Hcurl_Space>,  Shape_Hash> shape_Hcurl_;
    std::vector<Basis_Shape> basis_shapes_;

public:
    Hcurl_Space() = default;
    Hcurl_Space(int dim, int p);

    bool add_basis_shape(Basis_Shape g) override;
    FEM_Space * get_basis_space(Basis_Shape s) const override;
    const std::vector<Basis_Shape>& get_basis_shapes() const override;


    // Returns basis values at a point in the unit tetrahedron
    // For H1: Scalars. For HCurl: Vectors.
    //void get_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const;
    //void get_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const;
    
    // Return vector proxy of Exterior Derivative of basis of the corresponding form.
    //void get_ED_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const;
    //void get_ED_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const;

    Space get_function_space() const override {return Space::H_curl;};

protected:
    // these function should not be called at this abstraction level!
    // following function cannot be called from Hcurl_space instance, it must be called from Hcurl_Tetrahedron for example.
    // For H1: Scalars. For HCurl: Vectors.
    virtual void get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const { Logger::error("Hcurl_Space::get_basis_s should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }
    virtual void get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const { Logger::error("Hcurl_Space::get_basis_v should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }

    // Return vector proxy of Exterior Derivative of basis of the corresponding form.
    virtual void get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const { Logger::error("Hcurl_Space::get_ED_basis_s should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }
    virtual void get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const { Logger::error("Hcurl_Space::get_ED_basis_v should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }

    virtual void dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const { Logger::error("Hcurl_Space::dof_transformation should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }

    virtual void dof_signature(
        int entity_dim, 
        int entity_idx, 
        const std::vector<Ref_Coord>& coord_list, 
        Eigen::Ref<MatrixXd> kernel
    ) const { Logger::error("Hcurl_Space::get_dof_signature should not be called at this abstraction level (please call from, e.g., Hcurl_triangle)."); }
};



/**
 * Reference tetrahedron:
 * 
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0] 
 *
 */
class Hcurl_triangle : public Hcurl_Space 
{
    friend class Hcurl_Space;
private:
    Hcurl_triangle(int p = 1);
    ~Hcurl_triangle(){};

public:

    void get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const override;
    void get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const override;

    void get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const override;
    void get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const override;

    void dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const override;

    void dof_signature(
        int entity_dim, 
        int entity_idx, 
        const std::vector<Ref_Coord>& coord_list, 
        Eigen::Ref<MatrixXd> kernel
    ) const override;
};



/**
 * Reference tetrahedron:
 * 
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0]  a3=[0,0,1]    
 * 
 */
class Hcurl_tetrahedron : public Hcurl_Space {
    friend class Hcurl_Space;
private:
    //std::vector<std::vector<double>> basis_cache;

    Hcurl_tetrahedron(int p = 1);
    ~Hcurl_tetrahedron(){};

public:
    //void get_basis_v(const Integration_Point p, const Eigen::Ref<MatrixXd> basis) const override {
    //    const double * ss = basis.data();
    //};


    void get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const override;
    void get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const override;

    void get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const override;
    void get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const override;

    void dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const override;

    void dof_signature(
        int entity_dim, 
        int entity_idx, 
        const std::vector<Ref_Coord>& coord_list, 
        Eigen::Ref<MatrixXd> kernel
    ) const override;

};


}