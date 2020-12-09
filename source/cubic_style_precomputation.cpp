#include "../include/cubic_style_precomputation.h"
#include "../include/edges.h"
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void cubic_style_precomputation(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    cubic_style_data & data)
{
    //get the cotangent matrix
    igl::cotmatrix(V, F, data.L);
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(data.L, b, Aeq, false, data.data);
    //get the barycentric area matrix
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, data.Area);
    
    //compute K: same as the implementation from the deformation assignment
    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V, F, C);
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(F.rows() * F.cols() * F.cols() * 3 * 2);
    
    int i, j, idx;
    for (int f = 0; f < F.rows(); f ++) {
        
        for (int e = 0; e < F.cols(); e ++) {
            i = F(f, (e + 1) % 3);
            j = F(f, (e + 2) % 3);
            for (int k = 0; k < F.cols(); k ++) {
                idx = F(f, k);
                for (int beta = 0; beta < 3; beta ++) {
                    triplets.push_back(Eigen::Triplet<double>(i, 3 * idx + beta, C(f, e) * (V(i, beta) - V(j, beta)) / 3.0));
                    triplets.push_back(Eigen::Triplet<double>(j, 3 * idx + beta, -1.0 * C(f, e) * (V(i, beta) - V(j, beta)) / 3.0));
                }
            }
        }
    }
    data.K = Eigen::SparseMatrix<double>(V.rows(), V.rows() * 3);
    data.K.setFromTriplets(triplets.begin(), triplets.end());
    
    //get per vertex area-weighted normal
    igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, data.N);
    
    //get the vertex triangle adjacency list, used to find spokes and rims edges, VI is useless.
    std::vector<std::vector<int>> VTadjacency;
    std::vector<std::vector<int>> VI;
    igl::vertex_triangle_adjacency(V.rows(), F, VTadjacency, VI);
    
    std::vector<int> triangle_adjacency;
    Eigen::MatrixXi Fi;
    Eigen::MatrixXd Di, Wi;
    data.E.resize(V.rows());
    data.D.resize(V.rows());
    data.W.resize(V.rows());
    for (int i = 0; i < V.rows(); i ++) {
        triangle_adjacency = VTadjacency[i];
        
        //construct Ei. note that there could be edge shared by two triangles. need to remove the duplicate.
        //IMPORTANT: Removing duplicate edge results in same result as including duplicate edge, but removing duplicate edge takes longer time
        //first construct per-vertex face
        
        //Fi.resize(triangle_adjacency.size(), 3);
        //for (int j = 0; j < triangle_adjacency.size(); j ++) {
            //Fi.row(j) = F.row(triangle_adjacency[j]);
        //}
        //find unique edge by the method in assignment 1
        //data.E[i] = edges(Fi);
        

        data.E[i].resize(triangle_adjacency.size() * 3, 2);
        for (int j = 0; j < triangle_adjacency.size(); j++)
        {
            for (int k = 0; k < 3; k ++) {
                data.E[i](3 * j + k, 0) = F(triangle_adjacency[j], k);
                data.E[i](3 * j + k, 1) = F(triangle_adjacency[j], (k + 1) % 3);
            }
        }
        
        //construce Di and Wi based on Ei;
        Di.resize(data.E[i].rows(), 3);
        Wi.resize(data.E[i].rows(), data.E[i].rows());
        Wi.setZero();
        for (int e = 0; e < data.E[i].rows(); e ++) {
            Di.row(e) = V.row((data.E[i])(e, 1)) - V.row((data.E[i])(e, 0));
            Wi(e, e) = data.L.coeff((data.E[i])(e, 0), (data.E[i])(e, 1));
        }
        data.D[i] = Di.transpose();
        data.W[i] = Wi;
    }
    // initialize Rho, Z, U;
    data.Rho.resize(V.rows());
    //value of Rho is found in the paper
    data.Rho.setConstant(1e-4);
    //Z, U are initialized to 0, accodring to the paper.
    data.Z.resize(V.rows(), 3);
    data.Z.setZero();
    data.U.resize(V.rows(), 3);
    data.U.setZero();
}
