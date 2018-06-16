//
//  fusion.hpp
//  Fusion
//
//  Created by Karan on 4/9/18.
//  Copyright � 2018 Karan. All rights reserved.
//

#ifndef fusion_hpp
#define fusion_hpp

#include <stdio.h>
#include "lib/Eigen/Dense"
#include "ekf.hpp"
#include "geo_ned.hpp"
#include <iostream>

using namespace geodectic_converter;

//Input data format assumed as - latitude, longitude, vel, psi_dot, accel, alt

enum class DataPointType
{
    IMU, GPS
};

/**
 * @brief Data interface for getting fusion actionable data from raw sensor measurements
 */
class DataPoint
{
public:

    /**
     * @brief Default constructor
     */
    DataPoint(bool verbose = false)
    :_initialized(false), _first_data_point(true)
    {
        _dx = 0;
        _dy = 0;
        _mx = 0;
        _my = 0;
        _ds = 0;

        _RadiusEarth = 6378388.0; //m
        _arc = 2.0 * M_PI * (_RadiusEarth + 230)/360.0; // degree
       this->verbose = verbose;
       if(this->verbose) std::cout << "     DATAPOINT: ----- Initialized.....\r\n";
    }

    /**
     * @brief Retrieves raw sensor data and stores it in private variables
     *
     * @param timestamp Current timestamp for the sensor data
     * @param data_type Data type: Either GPS: which includes GPS+IMU data or IMU: which only includes the IMU data
     * @param raw_data Raw sensor data
     */
     // TODO remove data_type as it is not used anymore
    void set(const long long timestamp, const DataPointType data_type, const Eigen::VectorXd& raw_data)
    {
        if(this->verbose) std::cout << "        DATAPOINT: ----- In set\r\n";
        _raw_data.resize(raw_data.size());
        _timestamp = timestamp;
        _raw_data = raw_data;
        _initialized = true;

        if(_first_data_point && raw_data(0)!=0.0 && raw_data(1)!=0.0)
        {
            _dx = 0;
            _dy = 0;
            _mx = 0;
            _my = 0;
            _prev_lat = raw_data(0);
            _prev_long = raw_data(1);
            _arc = 2.0 * M_PI * (_RadiusEarth + raw_data(5))/360.0;
            _first_data_point = false;
        }
      else if(!_first_data_point)
      {
          _arc = 2.0 * M_PI * (_RadiusEarth + raw_data(5))/360.0;
          _dx = _arc * cos(raw_data(0) * M_PI/180.0) * (raw_data(1) - _prev_long);
          _dy = _arc * (raw_data(0) - _prev_lat);
          _ds = sqrt(_dx * _dx + _dy * _dy);

          if(_ds == 0.0)
          {
              _data_type = DataPointType::IMU;
          }
          else
          {
              _data_type = DataPointType::GPS;
          }

          _mx += _dx; // cumulative sum
          _my += _dy; // cumulative sum
        if(this->verbose) std::cout << " Mx, My: " << _mx << ", " << _my << std::endl;
          _prev_lat = raw_data(0);
          _prev_long = raw_data(1);
      }
    }

    /**
     * @brief Returns saved raw data for sensor fusion
     *
     * @return Sensor data measurements
     */
    Eigen::VectorXd get_state() const
    {
        Eigen::VectorXd state(6);

        // Convert raw data to fusion readable states
        double x = _mx;
        double y = _my;
        double vel = _raw_data(2);
        double psi = 0;
        double psi_dot = _raw_data(3);
        double a = _raw_data(4);

        state << x, y, psi, vel, psi_dot, a;

        return state;
    }

    /**
     * @brief Get raw sensor data
     *
     * @return Raw sensor data
     */
    Eigen::VectorXd get_raw_data() const
    {
        return _raw_data;
    }

    /**
     * @brief Get data type associated with the data at current timestep
     *
     * @return Data type: Either GPS or IMU
     */
    DataPointType get_data_point_type() const
    {
        return _data_type;
    }

    /**
     * @brief Get current timestamp
     *
     * @return Timestamp associated with current data
     */
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

    double _prev_lat;
    double _prev_long;
    double _curr_lat;
    double _curr_long;
    double _dx;
    double _dy;
    double _mx;
    double _my;
    double _ds;

    double _RadiusEarth;
    double _arc;
    bool verbose;
};


/**
 * @brief Class which uses EKF and DataPoint to run the filter
 */
class Fusion
{
private:
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
    // GPS offsets
    double _xOffset, _yOffset;
    bool verbose;
public:

    /**
     * @brief Constructor
     *
     * @param max_acceleration System parameter specifying maximum acceleration
     * @param max_turn_rate System parameter specifying maximum turn rate
     * @param max_yaw_accel System parameter specifying max yaw acceleration
     */
    Fusion(double max_acceleration, double max_turn_rate, double max_yaw_accel, double varGPS,
    double varSpeed, double varYaw, double varAcc, double xOffset, double yOffset, bool verbose);

    /**
     * @brief Updates the Q matrix
     *
     * @param dt Timestep
     */
    void const updateQ(double dt);

    /**
     * @brief Starts the filter
     *
     * @param data Current sensor data
     */
    void start(const DataPoint& data);

    /**
     * @brief Performs the prediction and update step for the filter based on the datatype of the input data
     *
     * @param data Current sensor data
     */
    void compute(const DataPoint& data);

    /**
     * @brief Looping function to either start or run the filter in a loop
     *
     * @param data Current sensor data
     */
    void process(const DataPoint& data);

    /**
     * @brief Returns the estimated state of the system
     *
     * @return Vector containing the estimated states
     */
    Eigen::VectorXd get_resulting_state() const;
};

#endif /* fusion_hpp */
