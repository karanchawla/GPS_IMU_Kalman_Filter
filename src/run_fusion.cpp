//
//  run_fusion.cpp
//  Main Class running GPS INS Fusion
//
// Created by Karan on 5/2/18.
// Copyright � 2018 Karan. All rights reserved.

#include "fusion.hpp"
#include "csv_reader.hpp"
#include "run_fusion.hpp"
#include <eigen3/Eigen/Dense>

/**
 * @brief Constructor for the main class
 *
 * @param max_acceleration maximum acceleration for the system
 * @param max_turn_rate maximum turn rate for the system
 * @param max_yaw_accel maximum yaw acceleration for the system
 * These parameters are used to set up the system noise matrices in the EKF
 */
GpsIns::GpsIns(bool verbose)
{
	try
    {
        params = getDefaultParams();
    }
    catch (const ifstream::failure& e)
    {
        cout << "Exception opening/reading parameter file";
     }

    _raw_data = Eigen::VectorXd(6);
    filter = new Fusion(params.maxAccel, params.maxTurnRate, params.maxYawAccel, params.varGPS,
                        params.varSpeed, params.varYaw, params.varAcc, params.xOff, params.yOff, verbose);

    _sensor_data = new DataPoint(verbose);
    _prev_enc_counter = 0;
    _prev_gps_counter = 0;
    _prev_imu_counter = 0;
    _imucounter = 0;
    _gpscounter = 0;
    _enccounter = 0;
    _dt = 0;
    _prev_time = clock();
    _cur_time = 0;
    this->verbose = verbose;
    if(this->verbose) std::cout << "	GPS-INS: finished Initializing Constructor" << "\n";
    if(this->verbose){
	    std::printf("EKF Params Initialized:\r\n	maxAccel = %.3f\r\n	maxTurnRate = %.3f\r\n	maxYawAccel = %.3f\r\n	varGPS = %.3f\r\n	varSpeed = %.3f\r\n	 varYaw = %.3f\r\n	varAcc = %.3f\r\n", params.maxAccel, params.maxTurnRate, params.maxYawAccel, params.varGPS, params.varSpeed, params.varYaw, params.varAcc);

    }
}

/**
 * @brief Default destructor
 */
GpsIns::~GpsIns()
{}

/**
 * @brief Reads in the GPS data
 */
void GpsIns::read_gps_data(double lat, double lon, double alt)
{
	if(this->verbose) std::cout << "	GPS-INS: In Read GPS" << "\n";
	//m.lock();
    std::lock_guard<std::mutex> lock(m);
	_raw_data(0) = lat;
	_raw_data(1) = lon;
	_raw_data(5) = alt;
	_gpscounter += 1;
	//m.unlock();
	if(this->verbose) std::cout << "	GPS-INS: Raw GPS data --- " << _raw_data << "\n";
	if(this->verbose) std::cout << "	GPS-INS: Exiting read GPS" << "\n";
}

/**
 * @brief Reads IMU values
 *
 * @param yaw_rate psi_dot
 * @param long_accel longitudinal acceleration of the robot
 * @param timestamp current timestamp
 *
 */
void GpsIns::read_imu_data(double yaw_rate, double long_accel)
{
	if(this->verbose) std::cout << "	GPS-INS: In Read IMU" << "\n";
	//m.lock();
    std::lock_guard<std::mutex> lock(m);
	this->_cur_time = clock();
	this->_dt = _cur_time - _prev_time;

	_raw_data(3) = yaw_rate;
	_raw_data(4) = long_accel;
	_imucounter += 1;

	// update previous time only if dt > 0
	if(_dt > 0.0) this->_prev_time = _cur_time;
	//m.unlock();
	if(this->verbose) std::cout << "	GPS-INS: IMU Data -- " << _raw_data << "\n";
	if(this->verbose) std::cout << "	GPS-INS: Exiting read IMU" << "\n";
}

/**
 * @brief Reads in the encoders - called in tasks
 *
 * @param psi yaw values computed from encoders
 * @param vel velocity computed from encoders
 */
void GpsIns::read_encoders(double vel)
{
    if(this->verbose) std::cout << "	GPS-INS: Read encoder" << "\n";
    //m.lock();
    std::lock_guard<std::mutex> lock(m);
    _raw_data(2) = vel;
    _enccounter += 1;
    //m.unlock();
    if(this->verbose) std::cout << "	GPS-INS: Exiting Read encoder" << "\n";
}

/**
 * @brief Sets up the data in DataPoint class to be used up by the filter
 * How do you ensure data in this object is the current and has been completely filled at this timestep?
 */
void GpsIns::set_data()
{
    if(this->verbose) std::cout << "	GPS-INS: Setting data" << "\n";

    bool flag = (_gpscounter > _prev_gps_counter);
    // iF gps got updated then we update the filter else we just predict
    if(flag)
        _data_type = DataPointType::GPS;
    else
        _data_type = DataPointType::IMU;

    if(this->verbose) std::cout << "	GPS-INS: delta Time used --- " << (float)_dt << "\n";
    _sensor_data->set(_dt, _data_type, _raw_data);
    if(this->verbose) std::cout << "	GPS-INS: Data set" << "\n";

}

/**
 * @brief Main loop for the filter
 */
void GpsIns::loop()
{
    //m.lock();
    std::lock_guard<std::mutex> lock(m);
    set_data();
    filter->process(*_sensor_data);
    _prev_gps_counter = _gpscounter;
    //m.unlock();
}

/**
 * @brief Returns the estimated state of the system
 *
 * @return Estimated state for the robot
 */
Eigen::VectorXd GpsIns::get_estimated_state() const
{
    return filter->get_resulting_state();
}

DataPoint GpsIns::get_sensor_data()
{
    return *_sensor_data;
}
