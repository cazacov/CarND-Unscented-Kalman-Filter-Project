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

  // State dimension
  n_x_ = 5;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;   // Half of gravitational acceleration is a reasonable estimation of max. possible acceleration one can reach with muscle force

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;   // That's about 57 degrees / sec2.
  
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

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma points count
  points_count_ = n_aug_ * 2 + 1;

  // Sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;

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
    return;
  }

  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;  //delta_t - expressed in seconds
  time_us_= meas_package.timestamp_;

  Prediction(delta_t);

  switch (meas_package.sensor_type_){
    case MeasurementPackage::LASER:
      if (use_laser_)
      {
        UpdateLidar(meas_package);
      }
      break;
    case MeasurementPackage::RADAR:
      if (use_radar_)
      {
        UpdateRadar(meas_package);
      }
      break;
    default:
      throw std::invalid_argument("Unsupported sensor type");
  }

  // print the output
  cout << "new x_ " << endl << x_ << endl;
  cout << "new P_ " << endl << P_ << endl << endl << endl;

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
        0,  0,  5 * 5,  0, 0,             // 5 m/s = 18 kmh is a reasonable estimation of bicycle speed
        0,  0,  0,  M_PI * M_PI, 0,       // PI is our max. uncertainty in initial bicycle direction
        0,  0,  0,  0,  1 * 1;            // 1 radian / second is turn speed estimation

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;

  // process noise covariance matrix
  Q_ = MatrixXd(2,2);
  Q_ << std_a_ * std_a_, 0,
        0,  std_yawdd_ * std_yawdd_;


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // generate sigma points
  MatrixXd Xsig = MatrixXd(n_aug_, points_count_);
  GenerateSigmaPoints(Xsig);

  // predict sigma points
  SigmaPointPrediction(Xsig, delta_t);

  // predict state and state covariance
  PredictMeanAndCovariance();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Map predicted sigma points to lidar measurement space

  // Lidar has 2 coordinates
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, points_count_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  for (int i = 0; i < points_count_; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);

    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  //calculate mean predicted measurement

  z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < points_count_; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }


  //measurement covariance matrix S

  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < points_count_; i++)
  {
    VectorXd t = Zsig.col(i) - z_pred;
    S += weights_(i) * t * t.transpose();
  }

  // Add lidar measurement noise
  S += R_laser_;

  UpdateState(Zsig, S, z_pred, meas_package.raw_measurements_, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // Map predicted sigma points to radar measurement space

  // Radar has 3 coordinates
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, points_count_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  for (int i = 0; i < points_count_; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double pv = Xsig_pred_(2,i);
    double ppsi = Xsig_pred_(3,i);

    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rho_dot = 0;

    if (fabs(rho) > 1E-6) // avoid division by zero
    {
      rho_dot = (px * cos(ppsi) * pv + py * sin(ppsi) * pv) / rho;
    }
    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_dot;
  }

  //calculate mean predicted measurement

  z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < points_count_; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S

  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < points_count_; i++)
  {
    VectorXd t = Zsig.col(i) - z_pred;
    S += weights_(i) * t * t.transpose();
  }

  // Add radar measurement noise
  S += R_radar_;

  UpdateState(Zsig, S, z_pred, meas_package.raw_measurements_, n_z);
}

void UKF::GenerateSigmaPoints(MatrixXd &Xsig)
{

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.block(0,0,n_x_,n_x_) = P_;
  P_aug.block(n_x_,n_x_,2,2) = Q_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  double scale = sqrt(lambda_ + n_aug_);

  MatrixXd t =  scale * A;

  Xsig.col(0) = x_aug;
  Xsig.block(0,1, n_aug_, n_aug_) = t.colwise() + x_aug;
  Xsig.block(0, 1 + n_aug_, n_aug_, n_aug_) = (-t).colwise() + x_aug;

}

void UKF::SigmaPointPrediction(MatrixXd &XSig, double delta_t) {

  MatrixXd predictions = MatrixXd(n_x_, points_count_);

  for (int i=0; i< points_count_; i++)
  {
    double v = XSig(2, i);
    double psi = XSig(3, i);
    double psi_dot = XSig(4, i);
    double n_a = XSig(5, i);
    double n_psi_dot = XSig(6, i);

    //predict sigma points
    VectorXd inc = VectorXd(n_x_);
    VectorXd inc2 = VectorXd(n_x_);

    //avoid division by zero
    if(fabs(psi_dot) < 1E-3)
    {
      inc(0) = v * cos(psi) * delta_t;
      inc(1) = v * sin(psi) * delta_t;
      inc(2) = 0;
      inc(3) = psi_dot * delta_t;
      inc(4) = 0;
    }
    else
    {
      inc(0) = v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi));
      inc(1) = v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi));
      inc(2) = 0;
      inc(3) = psi_dot * delta_t;
      inc(4) = 0;
    }

    inc2(0) = 0.5 * delta_t * delta_t * cos(psi) * n_a;
    inc2(1) = 0.5 * delta_t * delta_t * sin(psi) * n_a;
    inc2(2) = delta_t * n_a;
    inc2(3) = 0.5 * delta_t * delta_t * n_psi_dot;
    inc2(4) = delta_t * n_psi_dot;


    //write predicted sigma points into right column
    predictions.col(i) = XSig.col(i).head(n_x_) + inc + inc2;

  }

  Xsig_pred_ = predictions;

}

void UKF::PredictMeanAndCovariance() {

  //create vector for weights
  weights_ = VectorXd(points_count_);

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //set weights
  for (int i = 0; i < points_count_; i++)
  {
    if (i)
    {
      weights_(i) = 1.0 / 2.0 / (lambda_ + n_aug_);
    }
    else {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    }
  }

  //predict state mean
  x  = VectorXd::Zero(n_x_);
  for (int i = 0; i < points_count_; i++)
  {
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P  = MatrixXd::Zero(n_x_, n_x_);
  for (int i = 0; i < points_count_; i++)
  {
    VectorXd t =  Xsig_pred_.col(i) - x;
    t(3) = NormalizeAngle(t(3));

    MatrixXd tt = t.transpose();

    P = P + weights_(i) * t * tt;
  }

  // write result
  x_ = x;
  P_ = P;
}

void UKF::UpdateState(MatrixXd &z_sig, MatrixXd &s, VectorXd &z_pred, VectorXd &z, int n_z) {


  // Create matrix for cross correlation between sigma points in state space and measurement space

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  for (int i = 0; i < points_count_; i++)
  {
    //residual
    VectorXd z_diff = z_sig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //residual
  VectorXd z_diff = z - z_pred;

  //calculate Kalman gain K;
  MatrixXd kalman_gain = Tc * s.inverse();

  //update state mean and covariance matrix
  x_ = x_ + kalman_gain * z_diff;

  P_ = P_ - kalman_gain * s * kalman_gain.transpose();
}

/**
 * Maps input angle in range -PI..+PI
 */
double UKF::NormalizeAngle(double a)
{
  while (a > M_PI) a-= 2*M_PI;
  while (a < -M_PI) a+= 2*M_PI;
  return a;
}
