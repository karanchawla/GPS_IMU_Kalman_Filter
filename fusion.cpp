//
//  fusion.cpp
//  Fusion
//
//  Created by Karan on 4/9/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include "fusion.hpp"

Fusion::Fusion(double max_acceleration, double max_turn_rate, double max_yaw_accel)
: _initialized(false), _max_turn_rate(max_turn_rate), _max_acceleration(max_acceleration), _max_yaw_accel(max_yaw_accel)
{
    // Initialize initial uncertainity P0
    _P = Eigen::MatrixXd(_n, _n);
    _P << 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1000.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1000.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1000.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1000.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1000.0;
    
    // Intialize dt or use log time?
    double varGPS = 5;
    double varSpeed = 0.1;
    double varYaw = 0.1;
    double varAcc = 1.0;
    _R = Eigen::MatrixXd(5, 5); //Assuming 5 sources of measurement
    _R << pow(varGPS, 2), 0.0, 0.0, 0.0, 0.0,
          0.0, pow(varGPS, 2), 0.0, 0.0, 0.0,
          0.0, 0.0, pow(varSpeed, 2), 0.0, 0.0,
          0.0, 0.0, 0.0, pow(varYaw, 2), 0.0,
          0.0, 0.0, 0.0, 0.0, pow(varAcc, 2);
}

void const Fusion::updateQ(double dt)
{
    // Process Noise Covariance Matrix Q
    _Q = Eigen::MatrixXd(_n, _n);
    _sGPS = 0.5 * _max_acceleration * pow(dt, 2);
    _sVelocity = _max_acceleration * dt;
    _sCourse = _max_turn_rate * dt;
    _sYaw = _max_yaw_accel * dt;
    _sAccel = _max_acceleration;
    _Q << pow(_sGPS, 2), 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, pow(_sGPS, 2), 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, pow(_sCourse, 2), 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, pow(_sVelocity, 2), 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, pow(_sYaw, 2), 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, pow(_sAccel, 2);
    
    _KF.setQ(_Q);
}

void Fusion::start(const DataPoint &data)
{
    _timestamp = data.get_timestamp();
    Eigen::VectorXd state = data.get_state();
    _KF.start(_n, state, _P, _F, _Q);
    _initialized = true;
}

void Fusion::compute(const DataPoint &data)
{
    /*******************************************
     * Prediction Step
     - Assumes current velocity is the same for this dt
     *******************************************/
    
    // Assuming 1.e6 for timestamp;
    const double dt = (data.get_timestamp() - _timestamp)/ 1.e6;
    _timestamp = data.get_timestamp();
    
    this->updateQ(dt);
    _KF.predict();
    
    /*******************************************
     * Update Step
     - Updates appropriate matrices given a measurement
     - Assumes measurement is received either from GPS or IMU
     *******************************************/
    const Eigen::VectorXd z = data.get_raw_data();
    const Eigen::VectorXd state = _KF.get_resulting_state();
    
    Eigen::VectorXd Hx;
    Eigen::MatrixXd JH;
    
    if(data.get_data_point_type() == DataPointType::GPS)
    {
        Eigen::VectorXd s = data.get_state();
        Hx << s(0),
              s(1),
              s(2),
              s(3),
              s(4);
        
        JH <<  1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    }
    else if(data.get_data_point_type() == DataPointType::IMU)
    {
        Eigen::VectorXd s = data.get_state();
        Hx << s(0),
              s(1),
              s(2),
              s(3),
              s(4);
        
        JH <<  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    }
    
    
    _KF.update(z, Hx, JH, _R);
}

void Fusion::process(const DataPoint &data)
{
    _initialized ? this->compute(data) : this->start(data);
}

Eigen::VectorXd Fusion::get_resulting_state() const
{
    return _KF.get_resulting_state();
}
