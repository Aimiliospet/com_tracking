#pragma once
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <memory>

using namespace Eigen;

class DJE {

//std::unique_ptr<Matrix<double, 15, 12>> Lcf;
/// Error Covariance, Linearized state transition model, Identity matrix,
/// state uncertainty matrix

std::unique_ptr<Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> P, Af, If, Qf, Kf, dxf, s, R, temp, z;
std::unique_ptr<MatrixXd> DARE(std::unique_ptr<MatrixXd>& F, std::unique_ptr<MatrixXd>& H, std::unique_ptr<MatrixXd>& Q, std::unique_ptr<MatrixXd>& R, int maxiterations, double tolerance);
std::unique_ptr<MatrixXd> computeDiscreteDyn(const std::unique_ptr<MatrixXd>& x_, const std::unique_ptr<MatrixXd>& jointvelocities);

/*
std::unique_ptr<Matrix<double, 15, 15>> P, Af, Acf, If;
/// Linearized Measurement model
std::unique_ptr<Matrix<double, 6, 15>> Hf, Hvf;
std::unique_ptr<Matrix<double, 3, 15>> Hv;
/// State-Input Uncertainty matrix
std::unique_ptr<Matrix<double, 12, 12>> Qf;
/// Kalman Gain
std::unique_ptr<Matrix<double, 15, 6>> Kf;
std::unique_ptr<Matrix<double, 15, 3>> Kv;
/// Correction state vector
std::unique_ptr<Matrix<double, 15, 1>> dxf;
/// Update error covariance and Measurement noise
std::unique_ptr<Matrix<double, 6, 6>> s, R;
/// position, velocity , acc bias, gyro bias, bias corrected acc, bias corrected gyr, temp vectors
std::unique_ptr<Vector3d> r, v, omega, f, fhat, omegahat, temp, omega_old;
/// Innovation vectors

std::unique_ptr<Matrix<double, 6, 1>> z;
std::unique_ptr<Vector3d> zv;
*/


public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        std::unique_ptr<Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x;
        
        //std::unique_ptr<Matrix<double, size*3, 1>> x;
        bool firstrun;

        int size;
        std::unique_ptr<Eigen::Matrix3d> Rib;
        std::unique_ptr<Eigen::Affine3d> Tib;
        std::unique_ptr<Eigen::Quaterniond> qib;
        double tolerance;
        int maxiterations;
        double dt;
        DJE();

		void setnoj(int sizeofjoints)
		{
			size = sizeofjoints;
            //x->resize(size, 1); maybe resize not needed and it might work with dynamic allocation 
		}

        std::unique_ptr<MatrixXd> DARE(std::unique_ptr<MatrixXd>& F, std::unique_ptr<MatrixXd>& H, std::unique_ptr<MatrixXd>& Q, std::unique_ptr<MatrixXd>& R, int maxiterations, double tolerance);
        std::unique_ptr<MatrixXd> DJE::computeDiscreteDyn(const std::unique_ptr<MatrixXd>& x_, const std::unique_ptr<MatrixXd>& jointvelocities);

		void predict(const std::unique_ptr<MatrixXd>& jointvelocities);
		
		/** @fn void updateWithTwistRotation(Vector3d y,Quaterniond qy);
		 *  @brief realises the  update step of the Error State Kalman Filter (ESKF) with a base linear velocity measurement and orientation measurement
		 *  @param y 3D base velociy measurement in the world frame
		 * 	@param qy orientation of the base w.r.t the world frame in quaternion
		 */
		void update(const std::unique_ptr<MatrixXd>& jointpositions);
		/**
		 *  @fn void init()
		 *  @brief Initializes the Base Estimator
		 *  @details
		 *   Initializes:  State-Error Covariance  P, State x, Linearization Matrices for process and measurement models Acf, Lcf, Hf and rest class variables
		 */
		void init();


};