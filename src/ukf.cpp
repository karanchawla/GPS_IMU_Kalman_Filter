//
//  ukf.cpp
//  EKF
//
//  Created by Karan on 4/24/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include "ukf.hpp"

UKF::UKF(const int nin, Eigen::MatrixXd& process_noise, Eigen::MatrixXd& initial_state, Eigen::MatrixXd& intial_covar, const double alpha, double beta, const double kappa)
:_num_states(nin), _num_sigmas(1 + _num_states * 2)
{
    _Q = process_noise;
    _state = initial_state;
    _P = intial_covar;
    _beta = beta;
    _alpha = alpha;
    _kappa = kappa;
    _initialized = true;
    
    _lambda = std::pow(_alpha, 2) * (_num_states + _kappa) - _num_states;
    _covariance_weights = Eigen::VectorXd(_num_sigmas);
    _mean_weights = Eigen::VectorXd(_num_sigmas);
    
    _covariance_weights(0) = (_lambda / (_num_states + _lambda)) + (1 - std::pow(_alpha, 2) + _beta);
    _mean_weights(0) = (_lambda / (_num_states + _lambda));
    
    for(int i = 1; i < _num_sigmas; i++)
    {
        _covariance_weights(i) = 1 / (2 * (_num_states + _lambda));
        _mean_weights(i) = 1 / (2 * (_num_states + _lambda));
    }
    
    _sigmas = get_sigmas(); // num_states x num_sigmas
}

Eigen::MatrixXd UKF::get_sigmas()
{
    // Probably a bug here, matrix manipulation is mind fuck
    Eigen::MatrixXd ret = Eigen::MatrixXd(_num_sigmas, _num_states);
    
    Eigen::MatrixXd temp_mat = (_num_states + _lambda) * _P;
    auto spr_mat = temp_mat.sqrt();
    ret.block(1, _num_states, 0, 0) = _state; // block of size 1 row and _num_states column at 0,0
    for(int i = 0; i < _num_states; i++)
    {
        ret.block(1, _num_states, i + 1, 0) = _state + spr_mat.block(1, _num_states, i, 0);
        ret.block(1, _num_states, i + _num_sigmas + 1, 0) = _state - spr_mat.block(1, _num_states, i, 0);
    }
    
    return ret.transpose();
}

void UKF::update(std::vector<int> indices, Eigen::VectorXd& data, Eigen::MatrixXd& R)
{
    // Performs a measurement update
    // @param states: list of indices (zero-indexed) of which states were measured, that is, which are being updated
    
    // create y, sigmas of just the states that are being updated
    Eigen::MatrixXd y(indices.size(), _num_sigmas); // num_updates x num_sigmas
    for(auto i: indices)
    {
        y.block(i, 0, 1, _num_sigmas) = _sigmas.block(1, _num_sigmas, i, 0);
    }
    
    // create y_mean, the mean of just the states that are being updated
    Eigen::VectorXd y_mean(indices.size());
    for(auto row: indices)
    {
        y(row) = _state(row);
    }
    
    auto y_diff(y);
    auto x_diff(_sigmas);
    
    for(int i = 0; i < _num_sigmas; i++)
    {
        for(int j = 0; j < indices.size(); j++)
        {
            y_diff(j,i) -= y_mean(j);
        }
        
        for(int j = 0; j < _num_states; j++)
        {
            x_diff(j,i) -= _state(j);
        }
    }
    
    // Covariance of measurement
    Eigen::MatrixXd p_yy(indices.size(), indices.size());
    for(int i = 0; i < _num_sigmas; i++)
    {
        Eigen::VectorXd val = y_diff.block(i, 0, indices.size(), 1);
        p_yy.block(0, i, indices.size(), 1) = _covariance_weights(i) * val * val.transpose();
    }
    // Add measurement noise
    p_yy += R;
    
    Eigen::MatrixXd p_xy(_num_states, indices.size());
    for(int i = 0; i < _num_sigmas; i++)
    {
        Eigen::VectorXd valx = y_diff.block(i, 0, indices.size(), 1);
        Eigen::VectorXd valy = x_diff.block(i, 0, _num_states, 1);
        p_xy.block(0, i, _num_states, 1) = _covariance_weights(i) * valx * valy.transpose();
    }
    
    auto K = p_xy * p_yy.inverse();
    auto y_actual = data;
    _state = _state + K * (y - y_actual);
    _P = _P - K * (p_yy) * K.transpose();
    _sigmas = get_sigmas();
}


void UKF::predict(double timestep)
{
    // Performs prediction step
    
    Eigen::MatrixXd trans_sigmas = _sigmas.transpose(); // num_sigmas x num_states
    
    for(int row = 0; row < trans_sigmas.rows(); row++)
    {
        // Get each row the transposed sigmas matrix
        Eigen::VectorXd temp_row = Eigen::VectorXd(trans_sigmas.block(row, 0, 1, _num_states));
        trans_sigmas.block(row, 0, 1, _num_states) = iterate_system_dynamics(temp_row, timestep);
    }
    
    _sigmas = trans_sigmas.transpose(); // num_states x num_sigmas
    
    Eigen::VectorXd x_out = Eigen::VectorXd(_state);
    
    for(int i = 0; i < _num_states; i++)
    {
        // the mean of that variable is the sum of
        // the weighted values of that variable for each iterated sigma point
        for(int j = 0; j < _num_sigmas; j++)
        {
            x_out(i) += _mean_weights(i) * _sigmas(i,j);
        }
    }
    
    Eigen::MatrixXd p_out = Eigen::MatrixXd(_num_states, _num_states);
    // For each sigma point
    for(int i = 0; i < _num_sigmas; i++)
    {
        //take the distance from the mean
        //make it a covariance by multiplying by the transpose
        //weight it using the calculated weighting factor
        //and sum
        Eigen::VectorXd diff = trans_sigmas.block(1, _num_states, i, 0);
        p_out += _covariance_weights(i) * diff.transpose()* diff;
    }
    
    // Add process noise
    p_out += timestep * _Q;
    _P = p_out;
    _state = x_out;
}


Eigen::VectorXd UKF::iterate_system_dynamics(Eigen::VectorXd &x, double dt)
{
    Eigen::VectorXd ret = Eigen::VectorXd(x);
    ret(0) = x(0) + (x[3]/x[4]) * (sin(x(4) * dt + x(2)) - sin(x(2)));
    ret(1) = x(1) + (x(3)/x(4)) * (-cos(x(4) * dt + x(2))+ cos(x(2)));
    ret(2) = std::fmod(x(2) + x(4) * dt + M_PI, 2.0 * M_PI) - M_PI;
    ret(3) = x(3) + x(5) * dt;
    ret(4) = x(4) - x(6); // Constant Turn Rate
    ret(5) = x(5) - x(7); // Constant Acceleration
    
    return ret;
}
