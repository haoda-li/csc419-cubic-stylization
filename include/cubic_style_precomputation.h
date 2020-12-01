#ifndef CUBE_STYLE_PRECOMPUTATION_H
#define CUBE_STYLE_PRECOMPUTATION_H

#include <cubic_style_data.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

// Precompute data needed to efficiently conduct local-global iterations for
// cubic stylization.
//
// Inputs:
//   V  #V by 3 list of vertex positions
//   F  #F by 3 list of triangle indices into the rows of V
//   b  #b indices into V of handle vertices
//   data struct that may contain additional input constraints etc
// Outputs:
//   data struct that contains all the precomputed information for cubic stylization
void cubic_style_precomputation(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    cubic_style_data & data);

#endif
