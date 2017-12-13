cmake_minimum_required(VERSION 3.5)
project(Dependability)

set(CMAKE_CXX_STANDARD 14)

add_executable(Dependability main.cpp SystemStateVector.cpp SystemStateVector.h)