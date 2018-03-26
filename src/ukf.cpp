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
  std_a_ = 10;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 2 * n_x_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3.0 - n_x_;

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  if (!is_initialized_){
    Initialize(meas_package);
    is_initialized_ = true;
  }

  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

void UKF::Initialize(MeasurementPackage &meas_package)
{
  time_us_ = meas_package.timestamp_;

  double x = 0;
  double y = 0;

  switch (meas_package.sensor_type_)
  {
    case MeasurementPackage::LASER:
      x = meas_package.raw_measurements_[0]; // X
      y = meas_package.raw_measurements_[1]; // Y
      break;
    case MeasurementPackage::RADAR:
      {
        double rho = meas_package.raw_measurements_[0];
        double theta = meas_package.raw_measurements_[1];
        x = rho * cos(theta);
        y = rho * sin(theta);
      }
      break;
    default:
      throw std::invalid_argument("Unsupported sensor type");
  }

  //  Initialize state vector
  x_ = VectorXd(n_x_);
  x_ << x, y, 0, 0, 0;

  // state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1,  0,  0,  0,  0,
        0,  1,  0,  0,  0,
        0,  0,  std_a_ * std_a_, 0, 0,
        0,  0,  0,  M_PI * M_PI, 0,
        0,  0,  0,  0,  std_yawdd_ * std_yawdd_;

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  //measurement covariance matrix - radar
  R_radar_ << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;

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
}
