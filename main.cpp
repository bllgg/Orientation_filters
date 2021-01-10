#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

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
        float acc_x = (float)stof(row[0]);
        float acc_y = (float)stof(row[1]);
        float acc_z = (float)stof(row[2]);

        float gyr_x = (float)stof(row[3]);
        float gyr_y = (float)stof(row[4]);
        float gyr_z = (float)stof(row[5]);

        float mag_x = (float)stof(row[6]);
        float mag_y = (float)stof(row[7]);
        float mag_z = (float)stof(row[8]);

//        MadgwickAHRSupdateIMU(gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z);
        for (int i = 0; i < 10; i++) {
            sensor_fusion.MadgwickAHRSupdate(gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z, mag_x, mag_y, mag_z);
        }
        float roll = sensor_fusion.get_roll();
        float pitch = sensor_fusion.get_pitch();
        float yaw = sensor_fusion.get_yaw();
        float *current_quaternioun;
        current_quaternioun = sensor_fusion.get_q();
//        cout << roll << "," << pitch << "," << yaw <<endl;
        cout << current_quaternioun[0] << "," << current_quaternioun[1] << "," << current_quaternioun[2] << "," << current_quaternioun[3] << endl;
//        cout << current_quaternioun[0] << "," << current_quaternioun[1] << "," << current_quaternioun[2] << "," << current_quaternioun[3] << "\n" <<endl;

    }
}
