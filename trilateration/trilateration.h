#ifndef __TRILATERATION_H__
#define __TRILATERATION_H__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


bool lineartrilateration(Eigen::VectorXf& y, Eigen::VectorXf& w, Eigen::VectorXf& x, 
					Eigen::MatrixXf& Xi, Eigen::MatrixXf& C);

Eigen::VectorXf nonlinearTrilateration(Eigen::VectorXf& _x0, Eigen::VectorXf& _y, 
					Eigen::MatrixXf& _Xi, Eigen::VectorXf& _w);

// TODO: trilateration with RANSAC
void ransacTrilateration(void);


#endif // __TRILATERATION_H__