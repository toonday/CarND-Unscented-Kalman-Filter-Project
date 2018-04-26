#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  is_initialized_ = false;
  n_x_ = 5;											// init state dimension
  n_aug_ = 7;										// init augmented state dimension
  lambda_ = 3 - n_x_;								// init spreading parameter
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);		// init matrix with predicted sigma points as columns
  weights_ = VectorXd(2 * n_aug_ + 1);				// init vector for weights
  time_us_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {VectorXd, int} reference to vector, index position to normalize.
 */
void UKF::NormalizeAngle(VectorXd &vec, int index) {
    //angle normalization
    while (vec(index)> M_PI) vec(index)-=(2.*M_PI);
    while (vec(index)<-M_PI) vec(index)+=(2.*M_PI);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  // check if sensor type should be ignored
  if ((use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
      (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    
    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {
	  // init state vector
      x_ << 1, 1, 1, 1, 0.1;

      // init state covariance matrix
      P_ << 0.15,    0, 0, 0, 0,
               0, 0.15, 0, 0, 0,
               0,    0, 1, 0, 0,
               0,    0, 0, 1, 0,
               0,    0, 0, 0, 1;

      // init timestamp
      time_us_ = meas_package.timestamp_;

      if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
		// init state vector
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      }
      else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Convert radar from polar to cartesian coordinates
        float rho = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
		
		// init state vector
        x_(0) = rho * cos(phi);
        x_(1) = rho * sin(phi);
      }

      is_initialized_ = true;

      return;
    }
	
    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
	// compute elapsed time btw measurements
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(delta_t);

	
    /*****************************************************************************
    *  Update
    ****************************************************************************/
	// update UKF based on sensor type
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  /*****************************************************************************
  *  Generate Sigma Points
  ****************************************************************************/
  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  // set lambda for non-augmented sigma points
  lambda_ = 3 - n_x_;
  
  // calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  A = std::sqrt(lambda_+n_x_) * A;

  // set sigma points as columns of matrix Xsig
  MatrixXd x_part_1 = MatrixXd(n_x_, n_x_);
  MatrixXd x_part_2 = MatrixXd(n_x_, n_x_);
  for (int i=0; i<A.cols(); ++i) {
      x_part_1.col(i) = x_ + A.col(i);
      x_part_2.col(i) = x_ - A.col(i);
  }
  Xsig << x_, x_part_1, x_part_2;

  /*****************************************************************************
  *  Augment Sigma Points
  ****************************************************************************/
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // set lambda for augmented sigma points
  lambda_ = 3 - n_aug_;

  // create augmented mean state
  x_aug << x_, 0, 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.block(0,0,5,5) << P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  A_aug = std::sqrt(lambda_+n_aug_) * A_aug;

  // create augmented sigma points
  MatrixXd x_aug_part_1 = MatrixXd(n_aug_, n_aug_);
  MatrixXd x_aug_part_2 = MatrixXd(n_aug_, n_aug_);
  for (int i=0; i<A_aug.cols(); ++i) {
      x_aug_part_1.col(i) = x_aug + A_aug.col(i);
      x_aug_part_2.col(i) = x_aug - A_aug.col(i);
  }
  Xsig_aug << x_aug, x_aug_part_1, x_aug_part_2;

  /*****************************************************************************
  *  Predict Sigma Points
  ****************************************************************************/
  //predict sigma points
  for (int i=0; i<Xsig_pred_.cols(); ++i) {
      float phi_dot = Xsig_aug(4,i);
      float phi = Xsig_aug(3,i);
      float speed = Xsig_aug(2,i);
      float y_pos = Xsig_aug(1,i);
      float x_pos = Xsig_aug(0,i);
      
      if (phi_dot != 0) {
          Xsig_pred_.col(i) << x_pos + ((speed/phi_dot)*(sin(phi+(phi_dot*delta_t))-sin(phi))) + (0.5*delta_t*delta_t*cos(phi)*Xsig_aug(5,i)),
                              y_pos + ((speed/phi_dot)*(-cos(phi+(phi_dot*delta_t))+cos(phi))) + (0.5*delta_t*delta_t*sin(phi)*Xsig_aug(5,i)),
                              speed + 0 + (delta_t*Xsig_aug(5,i)),
                              phi + (phi_dot*delta_t) + (0.5*delta_t*delta_t*Xsig_aug(6,i)),
                              phi_dot + 0 + (delta_t*Xsig_aug(6,i));
      }
      else {
          Xsig_pred_.col(i) << x_pos + (speed*cos(phi)*delta_t) + (0.5*delta_t*delta_t*cos(phi)*Xsig_aug(5,i)),
                              y_pos + (speed*sin(phi)*delta_t) + (0.5*delta_t*delta_t*sin(phi)*Xsig_aug(5,i)),
                              speed + 0 + (delta_t*Xsig_aug(5,i)),
                              phi + 0 + (0.5*delta_t*delta_t*Xsig_aug(6,i)),
                              phi_dot + 0 + (delta_t*Xsig_aug(6,i));
      }
  }

  /*****************************************************************************
  *  Convert Predicted Sigma Points to Mean/Covariance
  ****************************************************************************/
  // set weights
  weights_.fill(0.0);
  weights_ = weights_.array() + (1 / (2*(lambda_+n_aug_)));
  weights_(0) = lambda_ / (lambda_+n_aug_);

  // predicted state mean
  x_.fill(0.0);
  x_ = Xsig_pred_ * weights_;

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
	NormalizeAngle(x_diff, 3);

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix irrespective of measurement type.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(VectorXd z,  bool isLidar, bool isRadar) {
  // set measurement dimension
  // lidar can measure p_x and p_y
  // radar can measure rho, phi, and rho_dot
  int n_z;
  if( isLidar ) {
	  n_z = 2;
  }
  else if( isRadar ) {
	  n_z = 3;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for(int i=0; i<Zsig.cols(); ++i) {
	  if( isLidar ) {
		  Zsig(0,i) = Xsig_pred_(0,i);
		  Zsig(1,i) = Xsig_pred_(1,i);
	  }
	  else if( isRadar ) {
		  Zsig(0,i) = sqrt(pow(Xsig_pred_(0,i),2) + pow(Xsig_pred_(1,i),2));
		  Zsig(1,i) = atan2(Xsig_pred_(1,i), Xsig_pred_(0,i));
		  Zsig(2,i) = ((Xsig_pred_(0,i)*cos(Xsig_pred_(3,i))*Xsig_pred_(2,i)) + (Xsig_pred_(1,i)*sin(Xsig_pred_(3,i))*Xsig_pred_(2,i))) / Zsig(0,i);
	  }
  }

  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  z_pred = Zsig * weights_;

  //calculate innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i<Zsig.cols(); ++i) {
      VectorXd z_diff = Zsig.col(i)-z_pred;
      S = S + (weights_(i) * (z_diff * z_diff.transpose()));
  }
  
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  if( isLidar ) {
	  R(0,0) = pow(std_laspx_, 2);
	  R(1,1) = pow(std_laspy_, 2);
  }
  else if( isRadar ) {
	  R(0,0) = pow(std_radr_, 2);
	  R(1,1) = pow(std_radphi_, 2);
	  R(2,2) = pow(std_radrd_, 2);
  }
  S = S + R;

  /*****************************************************************************
  *  UKF Update
  ****************************************************************************/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i<weights_.rows(); ++i) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
	  
	  if( isRadar ) {
		  NormalizeAngle(x_diff, 3);
		  NormalizeAngle(z_diff, 1);
	  }
      Tc = Tc + (weights_(i) * (x_diff * z_diff.transpose()));
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  
  if( isRadar ) {
	  NormalizeAngle(z_diff, 1);
  }

  //update state mean and covariance matrix
  x_ = x_ + (K * z_diff);
  P_ = P_ - (K * S * K.transpose());
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  Update(meas_package.raw_measurements_, true, false);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  Update(meas_package.raw_measurements_, false, true);
}
