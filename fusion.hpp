//
//  fusion.hpp
//  Fusion
//
//  Created by Karan on 4/9/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef fusion_hpp
#define fusion_hpp

#include <stdio.h>
#include "Eigen/Dense"
#include "ekf.hpp"
#include "geo_ned.hpp"


using namespace geodectic_converter;

enum class DataPointType
{
    IMU, GPS
};

class DataPoint
{
public:
    long long get_timestamp() const
    {
        return _timestamp;
    }
    
    Eigen::VectorXd get_state() const
    {
        Eigen::VectorXd state;
        // TO DO :
        // Convert raw data to fusion readable states
        
        return state;
    }
    
    Eigen::VectorXd get_raw_data() const
    {
        return _raw_data;
    }
    
    DataPointType get_data_point_type() const
    {
        return _data_type;
    }
    
    
private:
    long long _timestamp;
    Eigen::VectorXd _raw_data;
    DataPointType _data_type;
    
};

class Fusion
{
public:
    const int _n = 6;
    bool _initialized;
    long long timestamp;
    Eigen::MatrixXd _P;
    Eigen::MatrixXd _Q;
    Eigen::MatrixXd _F;
    Eigen::MatrixXd _R;
    EKF _KF;
    
    double _sGPS;
    double _sCourse;
    double _sVelocity;
    double _sYaw;
    double _sAccel;
    double _dt;
    long long _timestamp;
    
    double _max_turn_rate;
    double _max_acceleration;
    double _max_yaw_accel;
    
    
private:
    Fusion(double max_acceleration, double max_turn_rate, double max_yaw_accel);
    void const updateQ(double dt);
    void start(const DataPoint& data);
    void compute(const DataPoint& data);
    void process(const DataPoint& data);
    Eigen::VectorXd get_resulting_state() const;
    
    
};

#endif /* fusion_hpp */
