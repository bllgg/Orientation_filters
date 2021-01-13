//
// Created by bllgg on 1/10/21.
// Directly Inherited from Madgwick filter
//

#include "ahrs.h"
#include <math.h>
#include <stdlib.h>


/*==============================================================*/
/*                  Constructor and Destructor                  */
/*==============================================================*/
Madgwick::Madgwick(float beta_val){
    this->beta = beta_val;
}

Madgwick::~Madgwick(){

}

/*==============================================================*/
/*                       Public Functions                       */
/*==============================================================*/

void Madgwick::MadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float hx, hy;
    float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

    // Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
    if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
        MadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az);
        return;
    }

    gx = degrees_to_radians(gx);
    gy = degrees_to_radians(gy);
    gz = degrees_to_radians(gz);

    // Rate of change of quaternion from gyroscope
    qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

    // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
    if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

        // Normalise accelerometer measurement
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // Normalise magnetometer measurement
        recipNorm = invSqrt(mx * mx + my * my + mz * mz);
        mx *= recipNorm;
        my *= recipNorm;
        mz *= recipNorm;

        // Auxiliary variables to avoid repeated arithmetic
        _2q0mx = 2.0f * q0 * mx;
        _2q0my = 2.0f * q0 * my;
        _2q0mz = 2.0f * q0 * mz;
        _2q1mx = 2.0f * q1 * mx;
        _2q0 = 2.0f * q0;
        _2q1 = 2.0f * q1;
        _2q2 = 2.0f * q2;
        _2q3 = 2.0f * q3;
        _2q0q2 = 2.0f * q0 * q2;
        _2q2q3 = 2.0f * q2 * q3;
        q0q0 = q0 * q0;
        q0q1 = q0 * q1;
        q0q2 = q0 * q2;
        q0q3 = q0 * q3;
        q1q1 = q1 * q1;
        q1q2 = q1 * q2;
        q1q3 = q1 * q3;
        q2q2 = q2 * q2;
        q2q3 = q2 * q3;
        q3q3 = q3 * q3;

        // Reference direction of Earth's magnetic field
//        hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3;
//        hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3;
        hx = mx * (1 - 2 * (q0q0 + q1q1)) + my * 2 * (q1q2 + q0q3) + mz * 2 * (q1q3 + q0q2);
        hy = mx * 2 * (q1q2 - q0q3) + my * (1 - 2 * (q0q0 + q2q2)) + mz * 2 * (q2q3 + q0q1);
        _2bx = 2.0f * sqrt(hx * hx + hy * hy);
//        _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
        _2bz = 2.0f * (2 * mx * (q1q3 + q0q2) + 2 * my * (q2q3 + q0q1) + mz * (1 - 2 * (q0q0 + q3q3)));
        _4bx = 2.0f * _2bx;
        _4bz = 2.0f * _2bz;

        // Gradient decent algorithm corrective step
//        s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
//        s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
//        s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
//        s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s0 = -_2q2 * (2.0f * q1q3 - 2.0f * q0q2 - ax) + _2q1 * (2.0f * q0q1 + 2.0f * q2q3 - ay)                                                       - _2bz * q2                * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2               * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s1 =  _2q3 * (2.0f * q1q3 - 2.0f * q0q2 - ax) + _2q0 * (2.0f * q0q1 + 2.0f * q2q3 - ay) - 4.0f * q1 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3                * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + ( _2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s2 = -_2q0 * (2.0f * q1q3 - 2.0f * q0q2 - ax) + _2q3 * (2.0f * q0q1 + 2.0f * q2q3 - ay) - 4.0f * q2 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + ( _2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s3 =  _2q1 * (2.0f * q1q3 - 2.0f * q0q2 - ax) + _2q2 * (2.0f * q0q1 + 2.0f * q2q3 - ay)                                                       + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1               * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
        s0 *= -recipNorm;
        s1 *= -recipNorm;
        s2 *= -recipNorm;
        s3 *= -recipNorm;

        // Apply feedback step
        qDot1 -= beta * s0;
        qDot2 -= beta * s1;
        qDot3 -= beta * s2;
        qDot4 -= beta * s3;
    }

    // Integrate rate of change of quaternion to yield quaternion
    q0 += qDot1 * (1.0f / sampleFreq);
    q1 += qDot2 * (1.0f / sampleFreq);
    q2 += qDot3 * (1.0f / sampleFreq);
    q3 += qDot4 * (1.0f / sampleFreq);

    // Normalise quaternion
    recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;
}

void Madgwick::MadgwickAHRSupdateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

    gx = degrees_to_radians(gx);
    gy = degrees_to_radians(gy);
    gz = degrees_to_radians(gz);

    // Rate of change of quaternion from gyroscope
    qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

    // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
    if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

        // Normalise accelerometer measurement
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // Auxiliary variables to avoid repeated arithmetic
        _2q0 = 2.0f * q0;
        _2q1 = 2.0f * q1;
        _2q2 = 2.0f * q2;
        _2q3 = 2.0f * q3;
        _4q0 = 4.0f * q0;
        _4q1 = 4.0f * q1;
        _4q2 = 4.0f * q2;
        _8q1 = 8.0f * q1;
        _8q2 = 8.0f * q2;
        q0q0 = q0 * q0;
        q1q1 = q1 * q1;
        q2q2 = q2 * q2;
        q3q3 = q3 * q3;

        // Gradient decent algorithm corrective step
        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;

        // Apply feedback step
        qDot1 -= beta * s0;
        qDot2 -= beta * s1;
        qDot3 -= beta * s2;
        qDot4 -= beta * s3;
    }

    // Integrate rate of change of quaternion to yield quaternion
    q0 += qDot1 * (1.0f / sampleFreq);
    q1 += qDot2 * (1.0f / sampleFreq);
    q2 += qDot3 * (1.0f / sampleFreq);
    q3 += qDot4 * (1.0f / sampleFreq);

    // Normalise quaternion
    recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;
}

//void Madgwick::compute_orientation(float *q){
//    yaw  = radians_to_degrees(atan2(2*q[1]*q[2] + 2*q[0]*q[3], q[0]*q[0] + q[1]*q[1] - q[2]*q[2] -q[3]*q[3]));
//    pitch = radians_to_degrees(-1*asin(2*(q[1]*q[3] - q[0]*q[2])));
//    roll = radians_to_degrees(atan2(2*q[0]*q[1] + 2*q[2]*q[3], q[0]*q[0] + q[3]*q[3] - q[1]*q[1] - q[2]*q[2]));
//}
//
//float* Madgwick::quaternionMul(float *q_1, float *q_2){
//    float e_0, e_1, e_2, e_3;
//    e_0 = q_1[0]*q_2[0] - q_1[1]*q_2[1] - q_1[2]*q_2[2] - q_1[3]*q_2[3];
//    e_1 = q_1[0]*q_2[1] - q_1[1]*q_2[0] - q_1[2]*q_2[3] - q_1[3]*q_2[2];
//    e_2 = q_1[0]*q_2[2] - q_1[1]*q_2[3] - q_1[2]*q_2[0] - q_1[3]*q_2[1];
//    e_3 = q_1[0]*q_2[3] - q_1[1]*q_2[2] - q_1[2]*q_2[1] - q_1[3]*q_2[0];
//
//    float result_quat[4] = {e_0, e_1, e_2, e_3};
//    return result_quat;
//}
//
//// should clear pointer after using this function
//float* Madgwick::get_accel_jacobian(float *q){
//    float *jacob_mat = (float *)malloc(3 * 4 * sizeof(float));
//
//    jacob_mat[0 ] = -2*q[2];
//    jacob_mat[1 ] =  2*q[3];
//    jacob_mat[2 ] = -2*q[0];
//    jacob_mat[3 ] =  2*q[1];
//
//    jacob_mat[4 ] =  2*q[1];
//    jacob_mat[5 ] =  2*q[0];
//    jacob_mat[6 ] =  2*q[3];
//    jacob_mat[7 ] =  2*q[2];
//
//    jacob_mat[8 ] =       0;
//    jacob_mat[9 ] = -4*q[1];
//    jacob_mat[10] = -4*q[2];
//    jacob_mat[11] =       0;
//
//    return jacob_mat;
//}

// should clear pointer after using this function
//float* Madgwick::get_accel_function(float *q, float *a){
//    float *function = (float *)malloc(3 * sizeof(float));
//
//    function[0] = 2.0*(q[1]*q[3] - q[0]*q[2])       - a[1];
//    function[1] = 2.0*(q[0]*q[1] + q[2]*q[3])       - a[2];
//    function[2] = 2.0*(0.5 - q[1]*q[1] - q[2]*q[2]) - a[3];
//
//    return function;
//}
//
//void Madgwick::normalizeq(float *q){
//    float q_length = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
//    q[0] /= q_length;
//    q[1] /= q_length;
//    q[2] /= q_length;
//    q[3] /= q_length;
//}

//void Madgwick::update_roll_and_pitch(){
//
//}
//
//void Madgwick::get_mag_jacob(){
//
//}
//
//void Madgwick::get_mag_func(){
//
//}
//
//void Madgwick::get_rotation_mat(){
//
//}
//
//void Madgwick::update_roll_pitch_yaw(){
//
//}

/*==============================================================*/
/*                     Getters and Setters                      */
/*==============================================================*/

float Madgwick::get_roll(void){
    roll = radians_to_degrees(atan2(2 * this->q1 * this->q2 + 2 * this->q0 * this->q3, this->q0 * this->q0 + this->q1 * this->q1 - this->q2 * this->q2 - this->q3 * this->q3));
    return roll;
}

void Madgwick::set_roll(int roll){
    this->roll = roll;
}


float Madgwick::get_pitch(void){
    pitch = radians_to_degrees(-1 * asin(2 * (this->q1 * this->q3 - this->q0 * this->q2)));
    return pitch;
}

void Madgwick::set_pitch(int pitch){
    this->pitch = pitch;
}


float Madgwick::get_yaw(void){
    yaw = radians_to_degrees(atan2(2 * this->q1 * this->q2 + 2* this->q0 * this->q3, this->q0 * this->q0 + this->q1 * this->q1 - this->q2 * this->q2 - this->q3 * this->q3));
    return yaw;
}

void Madgwick::set_yaw(int yaw){
    this->yaw = yaw;
}


float Madgwick::get_beta(void){
    return beta;
}

void Madgwick::set_beta(float beta){
    this->beta = beta;
}


float* Madgwick::get_q(void){
    static float quaternion[4];
    quaternion[0] = q0;
    quaternion[1] = q1;
    quaternion[2] = q2;
    quaternion[3] = q3;

    return quaternion;
}

void Madgwick::set_q(float q_0, float q_1, float q_2, float q_3){
    q0 = q_0;
    q1 = q_1;
    q2 = q_2;
    q3 = q_3;
}

/*==============================================================*/
/*                      Private Functions                       */
/*==============================================================*/

float Madgwick::radians_to_degrees(float radian_val){
    float degree_val = radian_val*180/3.1415;
    return degree_val;
}

float Madgwick::degrees_to_radians(float degree_val) {
    float radian_val = degree_val*3.1415/180;
    return radian_val;
}

float Madgwick::invSqrt(float x) {
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long*)&y;
    i = 0x5f3759df - (i>>1);
    y = *(float *)&i;
    y = y * (1.5f - (halfx * y * y));
    return y;
}