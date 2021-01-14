#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "ahrs.h"
using namespace std;


int main() {
//    system("pwd");
    fstream file;
    file.open("/mnt/DAA45F54A45F31ED/Education/UOM_MY_WORK/FYP/Madgwick_Test/Datasets/IMU_Data_01_09_2021_16:10:33_roll.csv");
    if(!file.is_open())
    {
        printf("Error!");
        exit(1);
    }
    string line, word;
    vector<string> row;
    Madgwick sensor_fusion = Madgwick(0.5);

    while(file.good()){
        getline(file, line);
        // used for breaking words
        stringstream s(line);
        while (getline(s, word, ',')) {
            row.push_back(word);
        }
//        cout << row[0] << "," << row[1] << endl;
        double acc_x = (double)stod(row[0]);
        double acc_y = (double)stod(row[1]);
        double acc_z = (double)stod(row[2]);
        acc_x = (int)(round(acc_x * 1000)) / 1000.0;
        acc_y = (int)(round(acc_y * 1000)) / 1000.0;
        acc_z = (int)(round(acc_z * 1000)) / 1000.0;

        double gyr_x = (double)stof(row[3]);
        double gyr_y = (double)stof(row[4]);
        double gyr_z = (double)stof(row[5]);
        gyr_x = (int)(round(gyr_x * 1000)) / 1000.0;
        gyr_y = (int)(round(gyr_y * 1000)) / 1000.0;
        gyr_z = (int)(round(gyr_z * 1000)) / 1000.0;

        double mag_x = (double)stof(row[6]);
        double mag_y = (double)stof(row[7]);
        double mag_z = (double)stof(row[8]);
        mag_x = (int)(round(mag_x * 1000)) / 1000.0;
        mag_y = (int)(round(mag_y * 1000)) / 1000.0;
        mag_z = (int)(round(mag_z * 1000)) / 1000.0;

//        MadgwickAHRSupdateIMU(gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z);
        for (int i = 0; i < 1; i++) {
            sensor_fusion.MadgwickAHRSupdate(gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z, mag_x, mag_y, mag_z);
//            sensor_fusion.MadgwickAHRSupdateIMU(gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z);

        }
        double roll = sensor_fusion.get_roll();
        double pitch = sensor_fusion.get_pitch();
        double yaw = sensor_fusion.get_yaw();
        double *current_quaternioun;
        current_quaternioun = sensor_fusion.get_q();
//        cout << roll << "," << pitch << "," << yaw <<endl;
//        cout << roll << "," << pitch <<endl;

        cout << current_quaternioun[0] << "," << current_quaternioun[1] << "," << current_quaternioun[2] << "," << current_quaternioun[3] << endl;
//        cout << current_quaternioun[0] << "," << current_quaternioun[1] << "," << current_quaternioun[2] << "," << current_quaternioun[3] << "\n" <<endl;

    }
}
