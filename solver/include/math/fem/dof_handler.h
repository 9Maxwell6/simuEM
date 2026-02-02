#pragma once


class DoF_Handler{
    Mesh *mesh;
    const FiniteElementCollection *fec;

    size_t n_dof;
}