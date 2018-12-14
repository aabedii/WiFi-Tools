cmake_minimum_required(VERSION 3.5)
project(Transmission_time)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES_TT
    transmission-time.cpp)

set(SOURCE_FILES_E
    efficiency.cpp)

add_executable(Transmission_time ${SOURCE_FILES_TT})
add_executable(efficiency ${SOURCE_FILES_E})
