//
//  fusion.cpp
//  Fusion
//
//  Created by Karan on 4/9/18.
//  Copyright � 2018 Karan. All rights reserved.
//
#include <iostream>
#include "fusion.hpp"

Fusion::Fusion(double max_acceleration, double max_turn_rate, double max_yaw_accel, double varGPS,
    double varSpeed, double varYaw, double varAcc, double xOffset, double yOffset, bool verbose)
: _initialized(false), _max_turn_rate(max_turn_rate), _max_acceleration(max_acceleration), _max_yaw_accel(max_yaw_accel), _xOffset(xOffset),
  _yOffset(yOffset), _KF(verbose)
{
    // Initialize initial uncertainity P0
    _P = Eigen::MatrixXd(_n, _n);
    _P << 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1000.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1000.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1000.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1000.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1000.0;

    _R = Eigen::MatrixXd(5, 5); //Assuming 5 sources of measurement
    _R << pow(varGPS, 2), 0.0, 0.0, 0.0, 0.0,
          0.0, pow(varGPS, 2), 0.0, 0.0, 0.0,
          0.0, 0.0, pow(varSpeed, 2), 0.0, 0.0,
          0.0, 0.0, 0.0, pow(varYaw, 2), 0.0,
          0.0, 0.0, 0.0, 0.0, pow(varAcc, 2);

  this->verbose = verbose;
  if(verbose) std::cout << " =========================== FUSION:  Initializing --- " << "\r\n";
}

void const Fusion::updateQ(double dt)
{
  if(this->verbose) std::cout << " =========================== FUSION:  Updating Q --- " << "\r\n";

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
  if(this->verbose) std::cout << "    Fusion: ------ In start.....\r\n";
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
  if(this->verbose) std::cout << "    Fusion: ------ In compute.....\r\n";
    // Assuming 1.e6 for timestamp - confirm after running on the system
    const double dt = (data.get_timestamp())/ 1.e6;
    // const double dt = 0.1;
    if(this->verbose) std::cout << dt << "timestep in compute";
    _timestamp = data.get_timestamp();

    // Update Q
    this->updateQ(dt);
    // Update state and calculate jacobian
    _KF.updateJA(dt);
    // Prediction
    _KF.predict();

    /*******************************************
     * Update Step
     - Updates appropriate matrices given a measurement
     - Assumes measurement is received either from GPS or IMU
     *******************************************/
    Eigen::VectorXd zz = data.get_state();
    Eigen::VectorXd z;
    z.resize(5);
    z << zz(0), //east
         zz(1), //north
       zz(3), //vel
       zz(4), //yaw_rate
       zz(5); //accel

    const Eigen::VectorXd state = _KF.get_resulting_state();

    Eigen::VectorXd Hx;
    Eigen::MatrixXd JH;

  Hx.resize(5);
  JH.resize(5,6);
    
    // measurement function
    Hx << state(0) + _xOffset * cos(state(3)) - _yOffset * sin(state(3)),
          state(1) + _xOffset * sin(state(3)) + _yOffset * cos(state(3)),
          state(3),
          state(4),
          state(5);
          
    double j13 = - _xOffset * sin(state(3)) - _yOffset * cos(state(3));
    double j23 = _xOffset * cos(state(3)) - _yOffset * sin(state(3));
    if(data.get_data_point_type() == DataPointType::GPS)
    {
        JH <<  1.0, 0.0, j13, 0.0, 0.0, 0.0,
               0.0, 1.0, j23, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

        _KF.update(z, Hx, JH, _R);
    }
    else if(data.get_data_point_type() == DataPointType::IMU)
    {
        JH <<  0.0, 0.0, j13, 0.0, 0.0, 0.0,
               0.0, 0.0, j23, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
        _KF.update(z, Hx, JH, _R);
    }
}

void Fusion::process(const DataPoint &data)
{
  if(this->verbose) std::cout << "    Fusion: ------ In process.....\r\n";
    if(data.get_timestamp() > 0.0)
    {
        _initialized ? this->compute(data) : this->start(data);
    }
}

Eigen::VectorXd Fusion::get_resulting_state() const
{
    return _KF.get_resulting_state();
}
