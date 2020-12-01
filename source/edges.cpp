#include "edges.h"

Eigen::MatrixXi edges(const Eigen::MatrixXi &F)
{
    Eigen::MatrixXi E;
    // ADD YOUR CODE HERE
    
    Eigen::MatrixXi Edge_Checker = Eigen::MatrixXi::Zero(F.maxCoeff() + 1, F.maxCoeff() + 1);
    int counter = 0;
    int v_0, v_1, v_2;
    int e_0_start, e_0_end, e_1_start, e_1_end, e_2_start, e_2_end;
    
    for(int row = 0; row < F.rows(); row ++) {
        v_0 = F(row, 0);
        v_1 = F(row, 1);
        v_2 = F(row, 2);
        
        e_0_start = std::min(v_0, v_1);
        e_0_end = std::max(v_0, v_1);
        
        e_1_start = std::min(v_1, v_2);
        e_1_end = std::max(v_1, v_2);
        
        e_2_start = std::min(v_2, v_0);
        e_2_end = std::max(v_2, v_0);
        
        if (Edge_Checker(e_0_start, e_0_end) == 0) {
            Edge_Checker(e_0_start, e_0_end) += 1;
            counter += 1;
        }
        
        if (Edge_Checker(e_1_start, e_1_end) == 0) {
            Edge_Checker(e_1_start, e_1_end) += 1;
            counter += 1;
        }
        
        if (Edge_Checker(e_2_start, e_2_end) == 0) {
            Edge_Checker(e_2_start, e_2_end) += 1;
            counter += 1;
        }
    }
    E.resize(counter, 2);
    counter -= 1;
    
    for(int row = 0; row < F.rows(); row ++) {
        v_0 = F(row, 0);
        v_1 = F(row, 1);
        v_2 = F(row, 2);
        
        e_0_start = std::min(v_0, v_1);
        e_0_end = std::max(v_0, v_1);
        
        e_1_start = std::min(v_1, v_2);
        e_1_end = std::max(v_1, v_2);
        
        e_2_start = std::min(v_2, v_0);
        e_2_end = std::max(v_2, v_0);
        
        if (Edge_Checker(e_0_start, e_0_end) == 1) {
            Edge_Checker(e_0_start, e_0_end) -= 1;
            E(counter, 0) = e_0_start;
            E(counter, 1) = e_0_end;
            counter -= 1;
        }
        
        if (Edge_Checker(e_1_start, e_1_end) == 1) {
            Edge_Checker(e_1_start, e_1_end) -= 1;
            E(counter, 0) = e_1_start;
            E(counter, 1) = e_1_end;
            counter -= 1;
        }
        
        if (Edge_Checker(e_2_start, e_2_end) == 1) {
            Edge_Checker(e_2_start, e_2_end) -= 1;
            E(counter, 0) = e_2_start;
            E(counter, 1) = e_2_end;
            counter -= 1;
        }
    }
    return E;
}
