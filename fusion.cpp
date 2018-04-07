//
//  fusion.cpp
//  EKF
//
//  Created by Karan on 4/7/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include "fusion.hpp"

using namespace extended_kalman_filter;

EKF::EKF()
{
}

void EKF::start(const int nin, const VectorXd &xin, const MatrixXd &Pin, const MatrixXd &Fin, const MatrixXd &Qin)
{
    _num_states = nin;
    _I = MatrixXd::Identity(_num_states, _num_states);
    _state = xin;
    _P = Pin;
    _JA = Fin;
    _Q = Qin;
    
    return;
}
void EKF::setQ(const MatrixXd &Q_in)
{
    _Q = Q_in;
}

VectorXd EKF::get_resulting_state() const
{
    return _state;
}

void EKF::predict()
{
    _state(0) = _state(0) + (_state(3)/_state(4)) * (sin(_state(4) * _dt + _state(2)) - sin(_state(2)));
    _state(1) = _state(1) + (_state(3)/_state(4)) * (-cos(_state(4) * _dt + _state(2))+ cos(_state(2)));
    _state(2) = (_state(2) + _state(4) * _dt + M_PI) % (2.0 * M_PI) - M_PI;
    _state(3) = _state(3) + _state(5) * _dt;
    _state(4) = _state(4);
    _state(5) = _state(5);
    
    _JA = calculate_joacobian(_state, _dt);
}
