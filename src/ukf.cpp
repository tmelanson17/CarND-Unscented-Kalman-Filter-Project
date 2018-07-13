#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double normalizeAngle(double in){
  return max(-3.14159265, min(3.14159265, in));
}
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
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  n_x_ = 5;

  n_aug_ = 7;

  n_sigma_ = 2*n_aug_ + 1;
  
  lambda_ = 3 - n_x_;
  
  is_initialized_ = false;

  // Set vector for weights
  weights_ = VectorXd(n_sigma_);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_sigma_; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
  /**
  TODO:
    * Initialize the state x_ with the first measurement.
    * Create the covariance matrix.
    * Remember: you'll need to convert radar from polar to cartesian coordinates.
  */

  //state covariance matrix P
  P_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 5, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0.5;


  // first measurement
  cout << "UKF: " << endl;
  x_ = VectorXd(5); 
  x_ << 1, 1, 1, 1, 1; 

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    double rho = meas_package.raw_measurements_[0];
    double phi = meas_package.raw_measurements_[1];
    double rhodot = meas_package.raw_measurements_[2];
    x_ << rho*cos(phi), rho*sin(phi), rhodot, 0, 0; 
    
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
    */
    //set the state with the initial location and zero velocity
    x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

  }

  previous_timestamp_ = meas_package.timestamp_;
  // done initializing, no need to predict or update
  is_initialized_ = true;
  cout << "Done initializing" << endl;
  return;
  }

  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) return;
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) return;
  
  double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0 ;
  Prediction(delta_t);
  if (meas_package.sensor_type_  == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
  }
  // cout << "Finished update" << endl;
  previous_timestamp_ = meas_package.timestamp_;
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

  // Sigma estimation

  //create augmented mean state
  VectorXd x_aug(n_aug_);
  x_aug.setZero();
  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.setZero();
  MatrixXd Xsig_aug(n_aug_, n_sigma_);
  Xsig_aug.setZero();
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  Xsig_pred_.setZero();
  
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  // cout << "x_aug\n" << x_aug << endl;
  //create augmented covariance matrix
  P_aug.block(0, 0, n_x_, n_x_) = P_;
  P_aug.block(n_x_, n_x_, 2, 2) << std_a_*std_a_, 0,
                                 0, std_yawdd_*std_yawdd_;
  // cout << "P_aug\n" << P_aug << endl;
                        
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  MatrixXd cov_points = (sqrt(lambda_ + n_aug_) * A); // Positive end
  cov_points.colwise() += x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = cov_points;
  cov_points = (-sqrt(lambda_ + n_aug_) * A); // Negative end
  cov_points.colwise() += x_aug;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) = cov_points;
  //--

  // cout << "Xsig_aug\n" << Xsig_aug << endl;

  // Sigma Point Prediction 
  for (int i = 0; i < n_sigma_; i++){
    VectorXd x_i = Xsig_aug.col(i);
    double vk = x_i(2);
    double psik = x_i(3);
    double psidk = x_i(4);
    double nu_ak = x_i(5);
    double nu_yawddk = x_i(6);
    double epsilon = 0.001; // The error range for psidk == 0
    VectorXd noise(n_x_);
    noise << delta_t * delta_t / 2 * cos(psik) * nu_ak,
             0.5 * delta_t * delta_t * sin(psik) * nu_ak,
             delta_t * nu_ak,
             0.5 * delta_t * delta_t * nu_yawddk,
             delta_t * nu_yawddk;
    VectorXd delta_x(n_x_);
    if (abs(psidk) > epsilon){
        delta_x <<    vk / psidk * (sin(psik + psidk*delta_t) - sin(psik)),
                      vk / psidk * (-cos(psik + psidk*delta_t) + cos(psik)),
                      0,
                      psidk * delta_t,
                      0;
    } else {
        delta_x <<    vk * cos(psik) * delta_t,
                      vk * sin(psik) * delta_t,
                      0,
                      psidk * delta_t,
                      0;
    }
    if (i == 3){
    }
    Xsig_pred_.col(i) = x_i.head(n_x_) + delta_x + noise;
  }
  //--

  // cout << "Xsig_pred_\n" << Xsig_pred_ << endl;

  // Estimate new x, P values
  x_ = VectorXd(n_x_);
  x_.setZero();
  P_ = MatrixXd(n_x_, n_x_);
  P_.setZero();
  for (int i = 0; i < n_sigma_; i++){
    VectorXd xi = Xsig_pred_.col(i);
    //predict state mean
    // cout << "X vector " << xi.transpose() << " weight" << weights_(i) << endl;
    x_ += weights_(i) * xi;
  }
  for (int i = 0; i < n_sigma_; i++){
    //predict state covariance matrix
    P_ += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }
  // --

  // cout << "New x:\n" << x_ << endl;
  // cout << "New P:\n" << P_ << endl;
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

  int n_z = 2;
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd S(n_z, n_z);
  S.setZero();
  MatrixXd Zsig(n_z, Xsig_pred_.cols());
  Zsig.setZero();
  VectorXd z_pred(n_z);
  z_pred.setZero();

  //transform sigma points into measurement space
  for (int i=0; i<Xsig_pred_.cols(); i++){
      VectorXd xi = Xsig_pred_.col(i);
      double px = xi(0);
      double py = xi(1);
      VectorXd zi(n_z);
      zi << px, py; 
      Zsig.col(i) = zi;
  }

  
  //calculate mean predicted measurement
  for (int i=0; i<Zsig.cols(); i++){
      z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
  for (int i=0; i<Zsig.cols(); i++){
      VectorXd diff = (Zsig.col(i) - z_pred);
      diff(1) = normalizeAngle(diff(1));
      S += weights_(i) * diff * diff.transpose();
  }
  MatrixXd R(n_z, n_z);
  R.setZero();
  R(0,0) = std_laspx_ * std_laspx_;
  R(1,1) = std_laspy_ * std_laspy_;
  S += R;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  for (int i = 0; i < Xsig_pred_.cols(); i++){
      VectorXd diff = Zsig.col(i) - z_pred;
      VectorXd diffX = Xsig_pred_.col(i) - x_;
      diffX(3) = normalizeAngle(diffX(3));
      Tc += weights_(i) * (diffX) * (diff).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K*S*K.transpose();
 

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

  int n_z = 3;
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd S(n_z, n_z);
  S.setZero();
  MatrixXd Zsig(n_z, Xsig_pred_.cols());
  Zsig.setZero();
  VectorXd z_pred(n_z);
  z_pred.setZero();

  //transform sigma points into measurement space
  for (int i=0; i<Xsig_pred_.cols(); i++){
      VectorXd xi = Xsig_pred_.col(i);
      double px = xi(0);
      double py = xi(1);
      double v = xi(2);
      double psi = xi(3);
      VectorXd zi(n_z);
      zi << sqrt(px*px + py*py),
            atan2(py, px),
            (px*v*cos(psi) + py*v*sin(psi)) / sqrt(px*px + py*py);
      Zsig.col(i) = zi;
  }

  // cout << "Zsig : \n" << Zsig << endl;
  // cout << "z : \n" << z << endl;

  
  //calculate mean predicted measurement
  for (int i=0; i<Zsig.cols(); i++){
      z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
  for (int i=0; i<Zsig.cols(); i++){
      VectorXd diff = (Zsig.col(i) - z_pred);
      diff(1) = normalizeAngle(diff(1));
      // cout << "Diff: " << endl << diff.transpose() << endl;
      S += weights_(i) * diff * diff.transpose();
  }
  MatrixXd R(n_z, n_z);
  R.setZero();
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  S += R;

  // cout << "S matrix:\n" << S << endl;
  
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  for (int i = 0; i < Xsig_pred_.cols(); i++){
      VectorXd diff = Zsig.col(i) - z_pred;
      VectorXd diffX = Xsig_pred_.col(i) - x_;
      diff(1) = normalizeAngle(diff(1));
      diffX(3) = normalizeAngle(diffX(3));
      Tc += weights_(i) * (diffX) * (diff).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // cout << "Tc\n" << Tc <<endl;
  // cout << "K\n" << K << endl;

  //update state mean and covariance matrix
  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K*S*K.transpose();
 



}
