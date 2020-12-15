#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  R_laser_ << 0.0225, 0.0,
              0.0, 0.0225;

  R_radar_ << 0.09, 0.0, 0.0,
              0.0, 0.0009, 0.0,
              0.0, 0.0, 0.09;
  
  H_laser_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0;
  
  Hj_ << 1.0, 1.0, 0.0, 0.0,
         1.0, 1.0, 0.0, 0.0,
         1.0, 1.0, 1.0, 1.0; 
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  
  ekf_.F_ << 1.0, 0.0, 1.0, 0.0,
             0.0, 1.0, 0.0, 1.0,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0;

  ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1000.0, 0.0,
             0.0, 0.0, 0.0, 1000.0;
  

  ekf_.Q_ << 1.0, 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      long double rho = measurement_pack.raw_measurements_[0]; // range
      long double phi = measurement_pack.raw_measurements_[1]; // bearing
      long double rho_dot = measurement_pack.raw_measurements_[2]; // velocity
      
      long double x = rho * cos(phi);
      
      if ( x < 0.0001 ) {
        x = 0.0001;
      }
      
      long double y = rho * sin(phi);
      
      if ( y < 0.0001 ) {
        y = 0.0001;
      }      
      
      ekf_.x_(0) = rho * cos(phi);
      ekf_.x_(1) = rho * sin(phi);      
      ekf_.x_(2) = rho_dot * cos(phi);
      ekf_.x_(3) = rho_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;      
    }
    
  ekf_.F_ << 1.0, 0.0, 1.0, 0.0,
             0.0, 1.0, 0.0, 1.0,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0;

  ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1000.0, 0.0,
             0.0, 0.0, 0.0, 1000.0;  
    
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  
  // acceleration noise components 
  long double noise_ax = 9;
  long double noise_ay = 9;
  
  // time elapsed between the previous and current measurements
  long double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;  
  
  long double dt_2 = dt * dt;
  long double dt_3 = dt * dt * dt;
  long double dt_4 = dt * dt * dt * dt;
  
  //Modify the F matrix for time integration
	ekf_.F_ << 1, 0, dt, 0,
	          0, 1, 0, dt,
	          0, 0, 1, 0,
	          0, 0, 0, 1;

  if (dt > 0) {
    //process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_(0,0) = 0.250 * dt_4 * noise_ax;
    ekf_.Q_(0,1) = 0; 
    ekf_.Q_(0,2) = 0.50 * dt_3 * noise_ax;
    ekf_.Q_(0,3) = 0; 
    ekf_.Q_(1,0) = 0;    
    ekf_.Q_(1,1) = 0.250 * dt_4 * noise_ay;
    ekf_.Q_(1,2) = 0;   
    ekf_.Q_(1,3) = 0.50 * dt_3 * noise_ay;
    ekf_.Q_(2,0) = 0.50 * dt_3 * noise_ax;
    ekf_.Q_(2,1) = 0;
    ekf_.Q_(2,2) = dt_2 * noise_ax;
    ekf_.Q_(2,3) = 0;
    ekf_.Q_(3,0) = 0;
    ekf_.Q_(3,1) = 0.5 * dt_3 * noise_ay;
    ekf_.Q_(3,2) = 0;
    ekf_.Q_(3,3) = dt_2 * noise_ay;
  }
  
  ekf_.Predict();

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}