cmake_minimum_required(VERSION 3.6)
project(Transmission_time)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    transmission-time.cpp)

add_executable(Transmission_time ${SOURCE_FILES})