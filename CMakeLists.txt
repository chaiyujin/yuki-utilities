cmake_minimum_required(VERSION 2.6)

option(ENABLE_DEBUG "enable debug log" ON)
option(ENABLE_GUI "enable nanogui" ON)

# find necessary package
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${EIGEN3_INCLUDE_DIR})

# find gui
if (ENABLE_GUI)
    if (NOT EXISTS "${PROJECT_SOURCE_DIR}/ext/nanogui")
        message("Not found ext/nanogui git submodule, no gui related code will be compiled.")
    else()
        add_subdirectory(${PROJECT_SOURCE_DIR}/ext/nanogui/)
        add_definitions(-DWITH_NANOGUI)
    endif()
endif()

if (ENABLE_DEBUG)
else()
    add_definitions(-DNDEBUG)
    add_definitions(-O2)
endif()

add_definitions(-std=c++11)

add_executable(yuki
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/audio/io_wav.cpp
    ${PROJECT_SOURCE_DIR}/src/audio/features.cpp)
target_link_libraries(yuki Eigen3::Eigen fftw3)
