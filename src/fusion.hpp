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

//Input data format assumed as - latitude, longitude, psi, vel, psi_dot, accel, alt

enum class DataPointType
{
    IMU, GPS
};

class DataPoint
{
public:
    DataPoint()
    :_initialized(false), _first_data_point(true)
    {}
    
    void set(const long long timestamp, const DataPointType dta_type, const Eigen::VectorXd& raw_data)
    {
        _timestamp = timestamp;
        _data_type = _data_type;
        _raw_data = raw_data;
        _initialized = true;
        
        if(_first_data_point)
        {
            GPSDataConverter.intializeReference(raw_data(0), raw_data(1), raw_data(6));
            
            _first_data_point = false;
        }
        
        GPSDataConverter.geodetic2Ecef(raw_data(0), raw_data(1), raw_data(6), &_x, &_y, &_z);
        
        GPSDataConverter.ecef2Ned(_x, _y, _z, &_north, &_east, &_down);
    }
    
    Eigen::VectorXd get_state() const
    {
        Eigen::VectorXd state(6);
        
        // Convert raw data to fusion readable states
        if(_data_type == DataPointType::GPS)
        {
            double x = _north;
            double y = _east;
            double psi = _raw_data(2);
            double vel = _raw_data(3);
            double psi_dot = _raw_data(4);
            double a = _raw_data(5);
            
            state << x, y, psi, vel, psi_dot, a;
        }
        // If GPS data is not available
        else if(_data_type == DataPointType::IMU)
        {
            double x = 0.0;
            double y = 0.0;
            double psi = _raw_data(2);
            double vel = _raw_data(3);
            double psi_dot = _raw_data(4);
            double a = _raw_data(5);
            
            state << x, y, psi, vel, psi_dot, a;
        }
        
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
    
    long long get_timestamp() const
    {
        return _timestamp;
    }
    
private:
    double _x;
    double _y;
    double _z;
    double _north;
    double _east;
    double _down;
    GeodecticConverter GPSDataConverter;
    bool _initialized;
    bool _first_data_point;
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
    double _max_turn_rate;
    double _max_acceleration;
    double _max_yaw_accel;
    long long _timestamp;
    
private:
    Fusion(double max_acceleration, double max_turn_rate, double max_yaw_accel);
    void const updateQ(double dt);
    void start(const DataPoint& data);
    void compute(const DataPoint& data);
    void process(const DataPoint& data);
    Eigen::VectorXd get_resulting_state() const;
    
    
};

#endif /* fusion_hpp */
