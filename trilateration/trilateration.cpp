#include "trilateration.h"

#include <iostream>

using namespace Eigen;
using namespace std;

/*
 * y = Vector of N range measurements between beacon Xi and position X
 * w = Vector of associated range measurement variance
 * x = Position to estimate
 * Xi = Beacon position in matrix form, each row is a vector position
 * C = Covariance matrix of the estimated position X
 * */
bool linearTrilateration(VectorXf& y, VectorXf& w, VectorXf& x, 
					MatrixXf& Xi, MatrixXf& C) {

	unsigned int nRanges = y.size();
	unsigned int dim = x.size();
	unsigned int i = 0, j = 0;

	//cerr << "nRanges=" << nRanges;
	//cerr << ", dim=" << dim << ", x=" << x.transpose() << endl;

	if (((nRanges < 3) && (dim == 2)) || 
		((nRanges < 4) && (dim == 3))) {
		for (i=0;i<dim;i++) {
			x(i) = 1e12;
			for (j=0;j<dim;j++) {
				C(i,j) = 1e12;
			}
		}
		return false;
	}

	float si = 0.f;
	MatrixXf A = MatrixXf::Zero(nRanges, dim+1);
	MatrixXf invA, AAt, invG;
	MatrixXf W = MatrixXf::Zero(nRanges, nRanges);
	MatrixXf invW = MatrixXf::Zero(nRanges, nRanges);
	MatrixXf R = MatrixXf::Zero(nRanges, nRanges);
	MatrixXf G, T, C_temp;
	VectorXf b = VectorXf::Zero(nRanges);
	VectorXf d, x_temp;

	for (i=0;i<nRanges;i++) {
		si = 0.f;
		for (j=0;j<dim;j++) {
			A(i,j) = 2.f*Xi(i,j);
			si += pow(Xi(i,j),2.0f);
		}
		invW(i,i) = 1.f/w(i);
		if (!isfinite(invW(i,i)) || isnan(invW(i,i))) {
			invW(i,i) = 0.f;
		}
		W(i,i) = w(i);
		R(i,i) = -2.0f*y(i);
		A(i,dim) = -1.0f;
		b(i) = si - pow(y(i),2.0f);
	}
	//cerr << "w=" << w << endl;
	//cerr << "W=\n" << W << endl;
	//cerr << "A=\n" << A << endl;
	//cerr << "b=\n" << b << endl;
	G = A.transpose()*invW*A;
	d = A.transpose()*invW*b;
	AAt = A * A.transpose();
	invA = A.transpose()*AAt.inverse();
	invG = G.inverse();
	x_temp = invG*d;
	//cerr << "x_temp=\n" << x_temp << endl;
	for (i=0;i<dim;i++) {
		x(i) = x_temp(i);
		if (isnan(x(i)) || !isfinite(x(i))) {
			x(i) = 1e12;
			cerr << "trilateration: error: isnan or !isfinite\n";
			return false;
		}
	}

	/// Compute Covariance of the estimated pose
	T = invA*R;
	C_temp = T*W*T.transpose();
	//cerr << "C_temp=\n" << C_temp << endl;
	for (i=0;i<dim;i++) {
		for (j=0;j<dim;j++) {
			C(i,j) = C_temp(i,j);
			if (isnan(C(i,j)) || !isfinite(C(i,j))) {
				C(i,j) = 1000;
			}
		}
	}

	return true;
}


VectorXf nonlinearTrilateration(VectorXf& _x0, VectorXf& _y, 
					MatrixXf& _Xi, VectorXf& _w) {
	float min_e = 1.e-6;
	float diff_e = 1.e-6;
	float sum_e = 1.e6, sum_e_old = 1.e6;
	unsigned int max_iter = 400;
	unsigned int max_iter_follow = 10;
	unsigned int nDims = _x0.size();
	unsigned int nBeacons = _Xi.rows();
	unsigned int nRanges = _y.size();
	unsigned int iter = 0, iter_follow = 0;
	unsigned int i = 0, j = 0;
	VectorXf x = _x0;
	float d = 0.f, diff = 0.f;
	bool ok = true;

	if (nRanges != nBeacons)
		return x;

	VectorXf dx(nDims);
	VectorXf e = VectorXf::Zero(nBeacons);
	MatrixXf J = MatrixXf::Zero(nBeacons,nDims);
	MatrixXf Jt, gp, gpp, gppInv;
	VectorXf delta;

	do {
		sum_e = 0.f;
		for (i=0;i<nBeacons;i++) {
			d = 0.f;
			for (j=0;j<nDims;j++) {
				dx(j) = _Xi(i,j) - x(j);
				d += powf(dx(j),2.0f);
			}
			d = sqrtf(d);
			e(i) = d - _y(i);
			sum_e += e(i);
			d += (d == 0); // zero division
			for (j=0;j<nDims;j++)
				J(i,j) = -dx(j)/d;
		}
		Jt = J.transpose();
		gp = Jt*e;
		gpp = Jt*J;
		gppInv = gpp.inverse();
		delta = -gppInv*gp;
		x = x + delta;

		for (i=0;i<nDims;i++) {
			if (isnan(x(i)) || !isfinite(x(i))) {
				x(i) = 0;
				cerr << "nonlinearTrilateration: error, x isnan or !isfinite\n";
				ok = false;
			}
		}
		if (!ok)
			return x;

		sum_e = fabsf(sum_e);
		diff = fabsf(sum_e - sum_e_old);
		if (diff < diff_e)
			iter_follow++;
		else
			iter_follow = 0;
		sum_e_old = sum_e;
		iter++;
	} while (fabsf(sum_e)>min_e && iter<max_iter && iter_follow<max_iter_follow);

	return x;
}


void ransacTrilateration(void) {

}