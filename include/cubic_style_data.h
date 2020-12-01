#ifndef CUBE_STYLE_DATA_H
#define CUBE_STYLE_DATA_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>
struct cubic_style_data
{
    // Declare the members that can be precomputed
    // Write some brief documentation for each component you have defined
    //cotangent matrix L of the mesh
    Eigen::SparseMatrix<double> L;
    //store data used for global updates
    igl::min_quad_with_fixed_data<double> data;
    //barycentric area matrix Area
    Eigen::SparseMatrix<double> Area;
    //K: sparse matrix containing cotangents multiplied against differences across edges in the rest mesh.
    Eigen::SparseMatrix<double> K;
    //N: per vertex area-weighted normal;
    Eigen::MatrixXd N;
    //E: A list containing Ei for each vertex i. Ei contains all unique spoke and rim edges of vertex i, and each row of Ei contains the indices of end-points of a spoke edge or a rim edge. use Ei to construce Di~.
    std::vector<Eigen::MatrixXi> E;
    //D: A list containing Di for each vertex i. Di contains all unique spoke and rim edges of vertex i, and each row of Di represents a spoke edge or a rim edge.
    std::vector<Eigen::MatrixXd> D;
    //W: a list containing Wi for each vertex i. Wi is the per vertex cotangent matrix.
    std::vector<Eigen::MatrixXd> W;
    //Rho, Z, U: variables used in local update.
    Eigen::VectorXd Rho;
    Eigen::MatrixXd Z, U;
    
};

#endif
