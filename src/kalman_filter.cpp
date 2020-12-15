#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - (H_ * x_);
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd K = PHt * S.inverse();
  
  MatrixXd I = MatrixXd::Identity(4, 4);
  
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
    long double px = x_(0);
    long double py = x_(1);
    long double vx = x_(2);
    long double vy = x_(3);

    long double rho = sqrt(px*px + py*py);
    long double theta = atan2(py, px);
	long double rho_dot = (px*vx + py*vy) / rho;
  
    VectorXd h(3);
    h << rho, theta, rho_dot;
    VectorXd y = z - h;
 
      while (y(1)>M_PI) {
          y(1)-=2*M_PI;
      }

      while (y(1)<-M_PI) {
          y(1)+=2*M_PI;
      }
  
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;

    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
