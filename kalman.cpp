//
// Created by bllgg on 1/19/21.
//

#include <math.h>
#include <stdlib.h>
#include "kalman.h"

/*==============================================================*/
/*                  Constructor and Destructor                  */
/*==============================================================*/

Kalman::Kalman() = default;

Kalman::~Kalman() = default;

/*==============================================================*/
/*                       Public Functions                       */
/*==============================================================*/

void Kalman::KalmanFilter(double gx, double gy, double gz, double ax, double ay, double az){
    //prediction step
    theta_hat_roll = theta_t_minus_1_roll + delta_t * gx;
    theta_hat_pitch = theta_t_minus_1_pitch + delta_t * gy;

    p_hat_roll += gyro_noise;
    p_hat_pitch += gyro_noise;

    //measuring step
    z_k_roll = atan2(-ax, az); // roll angle around x axis
    z_k_pitch = atan2(-ay, sqrt(ax * ax + az * az)); // pitch angle around y axis

    //Estimation
    K_roll = p_hat_roll / (p_hat_roll + accel_noise);
    K_pitch = p_hat_pitch / (p_hat_pitch + accel_noise);

    theta_hat_roll  = theta_hat_roll + K_roll * (z_k_roll - theta_hat_roll);
    theta_hat_pitch = theta_hat_pitch + K_pitch * (z_k_pitch - theta_hat_pitch);

    p_hat_roll = (1 - K_roll) * p_hat_roll;
    p_hat_pitch = (1 - K_pitch) * p_hat_pitch;

}

/*==============================================================*/
/*                      Private Functions                       */
/*==============================================================*/

double Kalman::radians_to_degrees(double radian_val){
    double degree_val = radian_val * 180 / M_PI;
    return degree_val;
}

double Kalman::degrees_to_radians(double degree_val) {
    double radian_val = degree_val * M_PI / 180.0;
    return radian_val;
}

/*==============================================================*/
/*                     Getters and Setters                      */
/*==============================================================*/

double Kalman::get_roll() {
    roll = radians_to_degrees(theta_hat_roll);
    return roll;
}

double Kalman::get_pitch() {
    pitch = radians_to_degrees(theta_hat_pitch);
    return pitch;
}