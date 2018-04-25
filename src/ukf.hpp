//
//  ukf.hpp
//  EKF
//
//  Created by Karan on 4/24/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef ukf_hpp
#define ukf_hpp

#include <stdio.h>
#include "Eigen/Dense"
#include "Eigen/Core"
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <vector>

class UKF
{
public:
    UKF(const int nin, Eigen::MatrixXd& process_noise, Eigen::MatrixXd& initial_state, Eigen::MatrixXd& intial_covar, const double alpha, double beta, const double kappa);
    void update(std::vector<int> indices, Eigen::VectorXd& data, Eigen::MatrixXd& R);
    ~UKF(){};
    void predict(double timestep);
    Eigen::VectorXd& get_estimated_state() const;
    Eigen::VectorXd iterate_system_dynamics(Eigen::VectorXd &x, double dt);
private:
    const int _num_states;
    const int _num_sigmas;
    double _alpha;
    double _beta;
    double _kappa;
    double _lambda;
    bool _initialized;
    
    Eigen::VectorXd _covariance_weights;
    Eigen::VectorXd _mean_weights;
    Eigen::MatrixXd _sigmas;
    Eigen::MatrixXd _Q;
    Eigen::VectorXd _state;
    Eigen::MatrixXd _P;
    
    Eigen::MatrixXd get_sigmas();
};

#endif /* ukf_hpp */
