#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 1.5;

	/**
	 * DO NOT MODIFY measurement noise values below.
	 * These are provided by the sensor manufacturer.
	 */

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	 * End DO NOT MODIFY section for measurement noise values
	 */

	/**
	 * TODO: Complete the initialization. See ukf.h for other member properties.
	 * Hint: one or more values initialized above might be wildly off...
	 */

	// state vector
	n_x_ = 5;
	// augmented state vector
	n_aug_ = 7;
	// setup the augmented state vector
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	// lambda factor
	lambda_ = 3 - n_aug_;

	// setup the weights
	weights_ = VectorXd(2 * n_aug_ + 1);
	weights_.fill(0.5 / (lambda_ + n_aug_));
	weights_(0) = lambda_ / (lambda_ + n_aug_);

	// setup the noise covariance matrix
	R_lidar_ = MatrixXd(2, 2);
	R_lidar_.fill(0.0);
	R_lidar_(0, 0) = std_laspx_ * std_laspx_;
	R_lidar_(1, 1) = std_laspy_ * std_laspy_;

	R_radar_ = MatrixXd(3, 3);
	R_radar_.fill(0.0);
	R_radar_(0, 0) = std_radr_ * std_radr_;
	R_radar_(1, 1) = std_radphi_ * std_radphi_;
	R_radar_(2, 2) = std_radrd_ * std_radrd_;

	is_initialized_ = false;

	if (show_init_params_) {
		::std::cout << "------------------------- Init Parameters ----------------------" << std::endl;
		::std::cout << "std_a_ : " << std_a_ << ::std::endl;
		::std::cout << "std_yawdd_ : " << std_yawdd_ << ::std::endl;
		::std::cout << "R_radar_ : \n" << R_radar_ << ::std::endl;
		::std::cout << "R_lidar_ : \n" << R_lidar_ << ::std::endl;
		::std::cout << "----------------------------------------------------------------" << std::endl;
		show_init_params_ = false;
	}
}

void NormalizeAngle(double &angle) {
	while (angle > M_PI) angle -= 2. * M_PI;
	while (angle < -M_PI) angle += 2. * M_PI;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	 * TODO: Complete this function! Make sure you switch between lidar and radar
	 * measurements.
	 */

	if (!is_initialized_) {
		// Initialize the filter
		::std::cout << "----------------- Initializing the Kalman Filter ---------------" << std::endl;
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			std::cout << "Using the Radar Measurement Package !" << std::endl;
			double r = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rd = meas_package.raw_measurements_(2);
			double x = r * cos(phi);
			double y = r * sin(phi);
			double vx = rd * cos(phi);
			double vy = rd * sin(phi);
			double v = sqrt(vx * vx + vy * vy);
			x_ << x, y, 0.0, 0.0, 0.0;
			P_ = MatrixXd::Identity(n_x_, n_x_);
			P_(0, 0) = std_radr_ * std_radr_;
			P_(1, 1) = std_radr_ * std_radr_;
			P_(2, 2) = std_radrd_ * std_radrd_;
			P_(3, 3) = std_radphi_ * std_radphi_;
			P_(4, 4) = std_radphi_ * std_radphi_;

		} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			std::cout << "Using the Lidar Measurement Package !" << std::endl;
			double x = meas_package.raw_measurements_(0);
			double y = meas_package.raw_measurements_(1);
			x_ << x, y, 0.0, 0.0, 0.0;
			P_ = MatrixXd::Identity(n_x_, n_x_);
			P_(0, 0) = std_laspx_ * std_laspx_;
			P_(1, 1) = std_laspy_ * std_laspy_;
		}

		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		::std::cout << "Initializing Done !" << std::endl;
		::std::cout << "----------------------------------------------------------------" << std::endl;
		return;
	}

	// if already initialized, process the data normally
	::std::cout << "---------------- Processing the Measurement Data ---------------" << std::endl;
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
		// perform the predict step
		Prediction(dt);
		std::cout << "Prediction Done !" << std::endl;
		// perform the update step
		UpdateRadar(meas_package);
		std::cout << "Updated Radar Measurement !" << std::endl;

	} else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
		// perform the predict step
		Prediction(dt);
		std::cout << "Prediction Done !" << std::endl;
		// perform the update step
		UpdateLidar(meas_package);
		std::cout << "Updated Lidar Measurement !" << std::endl;
	}
	::std::cout << "----------------------------------------------------------------" << std::endl;
}

void UKF::Prediction(double delta_t) {
	/**
	 * TODO: Complete this function! Estimate the object's location.
	 * Modify the state vector, x_. Predict sigma points, the state,
	 * and the state covariance matrix.
	 */

	int c_z = 2 * n_aug_ + 1;
	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, c_z);
	Xsig_aug.fill(0.0);

	// get the augmented mean state
	x_aug.head(n_x_) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	// get the augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// save matrix
	MatrixXd L = P_aug.llt().matrixL();

	// get the augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int idx = 0; idx < n_aug_; ++idx) {
		Xsig_aug.col(idx + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(idx);
		Xsig_aug.col(idx + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(idx);
	}

	// predict sigma points
	for (int idx = 0; idx < c_z; ++idx) {
		double x = Xsig_aug(0, idx);
		double y = Xsig_aug(1, idx);
		double v = Xsig_aug(2, idx);
		double yaw = Xsig_aug(3, idx);
		double yaw_d = Xsig_aug(4, idx);
		double nu_a = Xsig_aug(5, idx);
		double nu_yawdd = Xsig_aug(6, idx);

		double x_p, y_p;

		if (fabs(yaw_d) > 0.001) {
			double yaw_p = yaw + yaw_d * delta_t;
			x_p = x + v / yaw_d * (sin(yaw_p) - sin(yaw));
			y_p = y + v / yaw_d * (cos(yaw) - cos(yaw_p));
		} else {
			x_p = x + delta_t * v * cos(yaw);
			y_p = y + delta_t * v * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yaw_d * delta_t;
		double yaw_d_p = yaw_d;

		x_p = x_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		y_p = y_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
		yaw_d_p = yaw_d_p + nu_yawdd * delta_t;

		Xsig_pred_(0, idx) = x_p;
		Xsig_pred_(1, idx) = y_p;
		Xsig_pred_(2, idx) = v_p;
		Xsig_pred_(3, idx) = yaw_p;
		Xsig_pred_(4, idx) = yaw_d_p;
	}

	// predict state
	VectorXd x_pred = VectorXd(n_x_);
	x_pred.fill(0.0);
	for (int idx = 0; idx < c_z; ++idx) {
		x_pred = x_pred + weights_(idx) * Xsig_pred_.col(idx);
	}
	x_ = x_pred;

	// predict covariance
	MatrixXd P_pred = MatrixXd(n_x_, n_x_);
	P_pred.fill(0.0);
	for (int idx = 0; idx < c_z; ++idx) {
		VectorXd x_diff = Xsig_pred_.col(idx) - x_;
		NormalizeAngle(x_diff(3));
		P_pred = P_pred + weights_(idx) * x_diff * x_diff.transpose();
	}
  P_ = P_pred;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	 * TODO: Complete this function! Use lidar data to update the belief
	 * about the object's position. Modify the state vector, x_, and
	 * covariance, P_.
	 * You can also calculate the lidar NIS, if desired.
	 */

	int n_z = 2;
	int c_z = 2 * n_aug_ + 1;
	VectorXd z = VectorXd(n_z);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);

	MatrixXd Zsig = MatrixXd(n_z, c_z);
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	// get the sigma points
	for (int idx = 0; idx < c_z; ++idx) {
		double x = Xsig_pred_(0, idx);
		double y = Xsig_pred_(1, idx);

		Zsig(0, idx) = x;
		Zsig(1, idx) = y;

		z_pred += weights_(idx) * Zsig.col(idx);
	}

	// get the matrix S
	for (int idx = 0; idx < c_z; ++idx) {
		VectorXd z_diff = Zsig.col(idx) - z_pred;
		NormalizeAngle(z_diff(1));
		S += weights_(idx) * z_diff * z_diff.transpose();
	}
	S += R_lidar_;

	// get the correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int idx = 0; idx < c_z; ++idx) {
		VectorXd z_diff = Zsig.col(idx) - z_pred;

		VectorXd x_diff = Xsig_pred_.col(idx) - x_;
		NormalizeAngle(x_diff(3));

		Tc += weights_(idx) * x_diff * z_diff.transpose();
	}

	//  get the Kalman gain K
	MatrixXd K = Tc * S.inverse();

	// update state vector and covariance matrix
	VectorXd z_diff = z - z_pred;
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	 * TODO: Complete this function! Use radar data to update the belief
	 * about the object's position. Modify the state vector, x_, and
	 * covariance, P_.
	 * You can also calculate the radar NIS, if desired.
	 */

	int n_z = 3;
	int c_z = 2 * n_aug_ + 1;
	VectorXd z = VectorXd(n_z);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
	z(2) = meas_package.raw_measurements_(2);

	MatrixXd Zsig = MatrixXd(n_z, c_z);
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	// get the sigma points
	for (int idx = 0; idx < c_z; ++idx) {
		double x = Xsig_pred_(0, idx);
		double y = Xsig_pred_(1, idx);
		double v = Xsig_pred_(2, idx);
		double yaw = Xsig_pred_(3, idx);
		double vx = cos(yaw) * v;
		double vy = sin(yaw) * v;

		Zsig(0, idx) = sqrt(x * x + y * y);
		Zsig(1, idx) = atan2(y, x);
		Zsig(2, idx) = (x * vx + y * vy) / sqrt(x * x + y * y);

		z_pred += weights_(idx) * Zsig.col(idx);
	}

	// get the matrix S
	for (int idx = 0; idx < c_z; ++idx) {
		VectorXd z_diff = Zsig.col(idx) - z_pred;
		NormalizeAngle(z_diff(1));
		S += weights_(idx) * z_diff * z_diff.transpose();
	}
	S += R_radar_;

	// get the correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int idx = 0; idx < c_z; ++idx) {
		VectorXd z_diff = Zsig.col(idx) - z_pred;
		NormalizeAngle(z_diff(1));

		VectorXd x_diff = Xsig_pred_.col(idx) - x_;
		NormalizeAngle(x_diff(3));

		Tc += weights_(idx) * x_diff * z_diff.transpose();
	}

	//  get the Kalman gain K
	MatrixXd K = Tc * S.inverse();

	// update state vector and covariance matrix
	VectorXd z_diff = z - z_pred;
	NormalizeAngle(z_diff(1));
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();
}