//
//  ekf.cpp
//  EKF
//
//  Created by Karan on 4/7/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include "ekf.hpp"

void EKF::start(const int nin, const Eigen::VectorXd &xin, const Eigen::MatrixXd &Pin, const Eigen::MatrixXd &Fin, const Eigen::MatrixXd &Qin)
{
    _num_states = nin;
    _I = Eigen::MatrixXd::Identity(_num_states, _num_states);
    _state = xin;
    _P = Pin;
    _JA = Fin;
    _Q = Qin;
    
    return;
}
void EKF::setQ(const Eigen::MatrixXd &Q_in)
{
    _Q = Q_in;
}
void EKF::updateJA(const double dt)
{
    // Updating state equations
    _state(0) = _state(0) + (_state(3)/_state(4)) * (sin(_state(4) * dt + _state(2)) - sin(_state(2)));
    _state(1) = _state(1) + (_state(3)/_state(4)) * (-cos(_state(4) * dt + _state(2)) + cos(_state(2)));
    _state(2) = std::fmod((_state(2) + _state(4) * dt + M_PI), (2.0 * M_PI)) - M_PI;
    _state(3) = _state(3) + _state(5) * dt;
    _state(4) = _state(4);
    _state(5) = _state(5);
    
    /* Calculate jacobian -
     linearizing the dynamics by 1st order Taylor series approximation */
    _JA = calculate_joacobian(_state, dt);
}

void EKF::predict()
{
    // Prediction step
    _P = _JA * _P * _JA.transpose() + _Q;
}

void EKF::update(const Eigen::VectorXd& Z, const Eigen::VectorXd& Hx, const Eigen::MatrixXd &JH, const Eigen::MatrixXd &R)
{
    // Tempoerary variable for storing this intermediate value
    _S = _JH * _P * _JH.transpose() + R;
    // Compute the Kalman gain
    _K = _P * _JH.transpose() * _S.inverse();
    // Update the estimate
    Eigen::VectorXd y = Z - Hx;
    _state = _state + _K * y;
    // Update the error covariance
    _P = (_I - _K * JH) * _P;
}

Eigen::VectorXd EKF::get_resulting_state() const
{
    return _state;
}
