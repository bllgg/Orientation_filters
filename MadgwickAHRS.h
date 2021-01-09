//
// Created by bllgg on 1/9/21.
//

#ifndef MADGWICK_TEST_MADGWICKAHRS_H
#define MADGWICK_TEST_MADGWICKAHRS_H

//----------------------------------------------------------------------------------------------------
// Variable declaration

extern volatile float beta;				// algorithm gain
extern volatile float q0, q1, q2, q3;	// quaternion of sensor frame relative to auxiliary frame

//---------------------------------------------------------------------------------------------------
// Function declarations

void MadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
void MadgwickAHRSupdateIMU(float gx, float gy, float gz, float ax, float ay, float az);

#endif //MADGWICK_TEST_MADGWICKAHRS_H
