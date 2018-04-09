//
//  ekf.hpp
//  EKF
//
//  Created by Karan on 4/7/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef ekf_hpp
#define ekf_hpp

#include <stdio.h>
#include "utils.hpp"
#include "Eigen/Dense"

class EKF
{
public:
    EKF(){};
    
    ~EKF(){};
    
    void start(const int nin, const Eigen::VectorXd& xin, const Eigen::MatrixXd& Pin, const Eigen::MatrixXd& Fin, const Eigen::MatrixXd& Qin);
    void setQ(const Eigen::MatrixXd& Q_in);
    Eigen::VectorXd get_resulting_state() const;
    void predict();
    void update(const Eigen::VectorXd& Z, const Eigen::VectorXd& Hx, const Eigen::MatrixXd &JH, const Eigen::MatrixXd &R);
private:
    bool _init;
    int _num_states;
    double _dt;
    double _dtGPS;
    
    Eigen::MatrixXd _P; // initial uncertainity
    Eigen::MatrixXd _Q; // process noise covariance
    Eigen::MatrixXd _JH; // measurment jacobian
    Eigen::MatrixXd _R; // measurement noise covariance
    Eigen::MatrixXd _I; // Identity matrix
    Eigen::MatrixXd _JA; // Jacobian state matrix
    Eigen::MatrixXd _S; // Matrix for storing intermediate step in update part
    Eigen::MatrixXd _K; // Kalman Gain
    Eigen::VectorXd _state; // 6 element state vector
};

#endif /* ekf_hpp */
