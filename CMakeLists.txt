cmake_minimum_required(VERSION 3.13)
project(Madgwick_Test)

set(CMAKE_CXX_STANDARD 14)

add_executable(Madgwick_Test main.cpp MadgwickAHRS.cpp MadgwickAHRS.h ahrs.cpp ahrs.h)