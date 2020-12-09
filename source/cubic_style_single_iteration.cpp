#include "../include/cubic_style_single_iteration.h"
#include "../include/edges.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/parallel_for.h>

void cubic_style_single_iteration(
	cubic_style_data& data,
	const Eigen::MatrixXd& bc,
	const double lambda,
	Eigen::MatrixXd& U)
{
	

	//local step first;
	//step 1: update R;
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3 * data.K.rows(), 3);
	
	// use igl's parallel to speed the local updates
	igl::parallel_for(data.K.rows(), [&data, &lambda, &R, &U](const int i) {
		//parameters for ADMM, found in the paper;
		double mu = 10.0;
		double t_incr = 2.0;
		double t_decr = 2.0;
		double e_abs = 1e-5;
		double e_rel = 1e-3;

		Eigen::MatrixXd M;
		Eigen::Matrix3d Ri;
		Eigen::MatrixXd D_tilda_i;

		D_tilda_i.resize(data.E[i].rows(), 3);
		for (int e = 0; e < data.E[i].rows(); e++) {
			D_tilda_i.row(e) = U.row((data.E[i])(e, 1)) - U.row((data.E[i])(e, 0));
		}
		//do local update until the norm of primal/dual residuals are small.
		while (true) {
			//step 1: update R
			//construct M
			M = data.D[i] * data.W[i] * D_tilda_i + data.Rho(i) * data.N.row(i).transpose() * (data.Z.row(i) - data.U.row(i));
			//SVD M
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			//make sure that the determinant of R is one
			Eigen::Matrix3d Omega;
			Omega.setIdentity(3, 3);
			Omega(2, 2) = (svd.matrixV() * (svd.matrixU().transpose())).determinant();
			Ri = svd.matrixV() * Omega * (svd.matrixU().transpose());

			//step 2: update z using shrinkage step
			//store previous z for calculating residual;
			Eigen::Vector3d z_prev = data.Z.row(i).transpose();
			Eigen::Vector3d x = Ri * data.N.row(i).transpose() + data.U.row(i).transpose();
			double k = lambda * data.Area.coeff(i, i) / data.Rho(i);
			for (int j = 0; j < 3; j++) {
				data.Z(i, j) = std::max((1.0 - k / std::abs(x(j))), 0.0) * x(j);
			}


			//step 3: get u~
			Eigen::Vector3d u_tilda;
			u_tilda = x - (data.Z.row(i).transpose());

			//stup 4: update u and rho;
			//calculate primal and dual residual;
			//formula of r and s are found in the paper, treat R as x and z as z; We need to flatten R such that we could get the correct A, B is a identity matrix;
			// I used row-wise flattened R, such that x.transpose = [R11, R12, R13, ..., R33].
			// r = Ax + Bz - C = -Rn + z
			double r = (data.Z.row(i).transpose() - Ri * data.N.row(i).transpose()).norm();
			// s = rho * AT * B * (z - z_prev) = rho * AT * (z - z_prev) since B is identity;
			//construct A first, not that we row-wise flattened R:
			//We could move A out of the while loop

			Eigen::MatrixXd A;
			A.resize(3, 9);
			A.setZero();
			A(0, 0) = -1.0 * data.N.row(i)(0);
			A(0, 1) = -1.0 * data.N.row(i)(1);
			A(0, 2) = -1.0 * data.N.row(i)(2);
			A(1, 3) = -1.0 * data.N.row(i)(0);
			A(1, 4) = -1.0 * data.N.row(i)(1);
			A(1, 5) = -1.0 * data.N.row(i)(2);
			A(2, 6) = -1.0 * data.N.row(i)(0);
			A(2, 7) = -1.0 * data.N.row(i)(1);
			A(2, 8) = -1.0 * data.N.row(i)(2);

			double s = (data.Rho(i) * A.transpose() * (data.Z.row(i).transpose() - z_prev)).norm();
			//updating rho and rescale U;
			if (r > mu * s) {
				data.Rho(i) = data.Rho(i) * t_incr;
				data.U.row(i) = u_tilda.transpose() / t_incr;
			}
			else if (s > mu * r) {
				data.Rho(i) = data.Rho(i) / t_decr;
				data.U.row(i) = u_tilda.transpose() * t_decr;
			}

			//check if we have reached the stopping criteria, i.e, if r and s are small:
			double e_primal = std::sqrt(6.0) * e_abs + e_rel * std::max((Ri * data.N.row(i).transpose()).norm(), (data.Z.row(i).transpose()).norm());
			//y = rho * u;
			double e_dual = std::sqrt(3.0) * e_abs + e_rel * (data.Rho(i) * A.transpose() * data.U.row(i).transpose()).norm();

			R.block(i * 3, 0, 3, 3) = Ri.transpose();
			if ((r < e_primal) && (s < e_dual)) {
				break;
			}

		}
		}
	);

	//global step then;
	Eigen::MatrixXd B = data.K * R;
	Eigen::MatrixXd Beq;
	igl::min_quad_with_fixed_solve(data.data, B, bc, Beq, U);
}
