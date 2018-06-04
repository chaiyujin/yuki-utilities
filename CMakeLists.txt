cmake_minimum_required(VERSION 2.6)

# find necessary package
find_package(Eigen3 3.3 REQUIRED)
include_directories(./src/)

# add definitions
add_definitions(-std=c++11 -O2)

add_executable(yuki
    ./src/main.cpp
    ./src/audio/io_wav.cpp)
target_link_libraries(yuki Eigen3::Eigen)