//
//  fusion.hpp
//  EKF
//
//  Created by Karan on 4/7/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef fusion_hpp
#define fusion_hpp

#include <stdio.h>
#include "utils.hpp"
#include "geo_ned.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace geodectic_converter;

namespace extended_kalman_filter
{
    class EKF
    {
    public:
        EKF();
        ~EKF();
        void start(const int nin, const VectorXd& xin, const MatrixXd& Pin, const MatrixXd& Fin, const MatrixXd& Qin);
        void initialize_states();
        void setQ(const MatrixXd &Q_in);
        VectorXd get_resulting_state() const;
        void predict();
    private:
        bool _init;
        int _num_states;
        double _dt;
        double _dtGPS;
        
        MatrixXd _P; // initial uncertainity
        MatrixXd _Q; // process noise covariance
        MatrixXd _JH; // measurment jacobian
        MatrixXd _R; // measurement noise covariance
        MatrixXd _I; // Identity matrix
        MatrixXd _JA; // Jacobian state matrix
        Eigen::VectorXd _state; // 6 element state vector
    };
};

#endif /* fusion_hpp */
