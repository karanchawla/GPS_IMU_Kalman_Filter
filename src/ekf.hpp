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
#include "lib/Eigen/Dense"

/**
 * @brief EKF base class implementing generic Extended Kalman Filter
 */
class EKF
{
public:

    /**
     * @brief Default constructor
     */
    EKF(bool verbose = false){this->verbose = verbose;};

    /**
     * @brief Default destructor
     */
    ~EKF(){};

    /**
     * @brief Function to be called the first time EKF is initialized
     *
     * @param nin Number of states
     * @param xin States
     * @param Pin Initial covariance matrix
     * @param Fin Initial Jacobian of the state
     * @param Qin Initial Q matrix
     */
    void start(const int nin, const Eigen::VectorXd& xin, const Eigen::MatrixXd& Pin, const Eigen::MatrixXd& Fin, const Eigen::MatrixXd& Qin);

    /**
     * @brief set Q value for the filter
     *
     * @param Q_in Input Q matrix
     */
    void setQ(const Eigen::MatrixXd& Q_in);

    /**
     * @brief Returns estimated state
     *
     * @return State of the system
     */
    Eigen::VectorXd get_resulting_state() const;

    /**
     * @brief Integrates system variables to predict the system state
     *
     * @param dt Time interval over which the integration takes place. Usually the difference between the previous and
     * current time step
     */
    void updateJA(const double dt);

    /**
     * @brief Updates the state covairance matrix and adds process noise
     */
    void predict();

    /**
     * @brief Runs the correction/update step of the filter
     *
     * @param Z Measurements for the current time step
     * @param Hx Measurement model
     * @param JH Jacobian of the measurment model
     * @param R Measurement noise
     */
    void update(const Eigen::VectorXd& Z, const Eigen::VectorXd& Hx, const Eigen::MatrixXd &JH, const Eigen::MatrixXd &R);
private:
    // Flag to indicate if the filter has started
    bool _init;
    bool verbose;
    int _num_states; // Number of states in the EKF

    Eigen::MatrixXd _P; // initial covaraince/uncertainity in states
    Eigen::MatrixXd _Q; // process noise covariance
    Eigen::MatrixXd _JH; // measurment jacobian
    Eigen::MatrixXd _R; // measurement noise covariance
    Eigen::MatrixXd _I; // Identity matrix
    Eigen::MatrixXd _JA; // Jacobian state matrix
    Eigen::MatrixXd _S; // Matrix for storing intermediate step in update part
    Eigen::MatrixXd _K; // Kalman Gain
    Eigen::VectorXd _state; // State - x y heading velocity yaw_rat long_acceleration
};

#endif /* ekf_hpp */
