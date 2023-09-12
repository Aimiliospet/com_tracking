#pragma once
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <memory>

using namespace Eigen;

class DVE {

//std::unique_ptr<Matrix<double, 15, 12>> Lcf;
/// Error Covariance, Linearized state transition model, Identity matrix,
/// state uncertainty matrix
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


/**
 * @brief computes the state transition matrix for linearized error state dynamics
 *
 */
std::unique_ptr<Matrix<double, 15, 15>> computeTrans(const std::unique_ptr<Matrix<double, 15, 1>>& x_,
													 const std::unique_ptr<Eigen::Matrix3d>& Rib_,
													 const std::unique_ptr<Vector3d>& omega_,
													 const std::unique_ptr<Vector3d>& f_);

/**
 * @brief performs euler (first-order) discretization to the nonlinear state-space dynamics
 *
 */
void euler(const std::unique_ptr<Vector3d>& omega_,const std::unique_ptr<Vector3d>& f_);

/**
 * @brief computes the discrete-time nonlinear state-space dynamics
 *
 */
std::unique_ptr<Matrix<double, 15, 1>> computeDiscreteDyn(const std::unique_ptr<Matrix<double, 15, 1>>& x_,
														  const std::unique_ptr<Eigen::Matrix3d>& Rib_,
														  const std::unique_ptr<Vector3d>& omega_,
														  const std::unique_ptr<Vector3d>& f_);

public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        std::unique_ptr<Matrix<double, 15, 1>> x;
        bool firstrun;
        std::unique_ptr<Vector3d> g, bgyr, bacc, gyro, acc, vel, pos, angle;
        double acc_qx, acc_qy, acc_qz, gyr_qx, gyr_qy, gyr_qz, gyrb_qx, gyrb_qy, gyrb_qz,
				accb_qx, accb_qy, accb_qz, odom_px, odom_py, odom_pz, odom_ax, odom_ay, odom_az,
				vel_px, vel_py, vel_pz, leg_odom_px, leg_odom_py, leg_odom_pz, leg_odom_ax,
				leg_odom_ay, leg_odom_az;

		double gyroX, gyroY, gyroZ, angleX, angleY, angleZ, bias_gx, bias_gy, bias_gz,
				bias_ax, bias_ay, bias_az, ghat;

		double accX, accY, accZ, velX, velY, velZ, rX, rY, rZ;
        std::unique_ptr<Eigen::Matrix3d> Rib;
        std::unique_ptr<Eigen::Affine3d> Tib;
        std::unique_ptr<Eigen::Quaterniond> qib;
        double dt;
        DVE();
        /** @fn void setdt(double dtt)
		 *  @brief sets the discretization of the Error State Kalman Filter (ESKF)
		 *  @param dtt sampling time in seconds
		 */
		void updateVars();


        /** @fn void predict(Vector3d omega_, Vector3d f_);
		 *  @brief realises the predict step of the Error State Kalman Filter (ESKF)
		 *  @param omega_ angular velocity of the base in the base frame
		 *  @param f_ linear acceleration of the base in the base frame
		 */
		void predict();

		 
		void update();
		/**
		 *  @fn void init()
		 *  @brief Initializes the Base Estimator
		 *  @details
		 *   Initializes:  State-Error Covariance  P, State x, Linearization Matrices for process and measurement models Acf, Lcf, Hf and rest class variables
		 */
		void init();

};