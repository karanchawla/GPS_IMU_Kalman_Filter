//
//  utils.cpp
//  EKF
//
//  Created by Karan on 4/9/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include "utils.hpp"

Eigen::MatrixXd calculate_joacobian(const Eigen::VectorXd& v, const double dt)
{
    // Assumes Jacobian is 6 x 6
    Eigen::MatrixXd JA = Eigen::MatrixXd::Zero(6,6);
    
    // Assumes the size of input vector is 6
    const double psi = v(2);
    const double velocity = v(3);
    const double psi_dot = v(4);
    
    const double THRESHOLD = 0.001;
    
    if (psi_dot > THRESHOLD) //Avoid dividing by zero
    {
        const double turn_radius = (velocity/psi_dot);
        const double psi_dot_inverse = 1/psi_dot;
        const double pdotp = dt * psi_dot + psi;
        
        const double r13 = turn_radius * (-cos(dt * psi_dot) + psi);
        const double r14 = psi_dot_inverse * (-sin(psi) + sin(pdotp));
        const double r15 = dt * turn_radius * cos(pdotp) - (turn_radius/psi_dot) * (-sin(psi) + sin(pdotp));
        
        const double r23 = turn_radius * (-sin(psi) + sin(pdotp));
        const double r24 = psi_dot_inverse * (cos(psi) - cos(pdotp));
        const double r25 = dt * turn_radius * sin(pdotp) - (turn_radius/psi_dot) * (cos(psi) - cos(pdotp));
        
        JA <<   1.0, 0.0, r13, r14, r15, 0.0,
        0.0, 1.0, r23, r24, r25, 0.0,
        0.0, 0.0, 1.0, 0.0, dt,  0.0,
        0.0, 0.0, 0.0, 1.0, 0.0,  dt,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ;
    }
    
    return JA;
}
