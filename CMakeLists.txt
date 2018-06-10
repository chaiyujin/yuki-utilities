cmake_minimum_required(VERSION 3.0)

option(ENABLE_DEBUG "enable debug log" OFF)
option(ENABLE_GUI "enable nanogui" OFF)

# find necessary package
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${EIGEN3_INCLUDE_DIR})

# find opencv
find_package(OpenCV REQUIRED)

# find openmp
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# add ffmpeg
if (APPLE)
    message("ffmpeg installed by homebrew")
    include_directories(/usr/local/Cellar/ffmpeg/4.0/include/)
    link_directories(/usr/local/Cellar/ffmpeg/4.0/lib/)

    set(FFMPEG_LIB avcodec avdevice avfilter avformat avresample avutil postproc swresample swscale)
endif()

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
    message("debug")
    add_definitions(-D_YUKI_DEBUG_)
else()
    message("release")
    add_definitions(-O2)
endif()

add_definitions(-std=c++11)

add_executable(yuki
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/math/mathutils.cpp
    ${PROJECT_SOURCE_DIR}/src/audio/io_wav.cpp
    ${PROJECT_SOURCE_DIR}/src/audio/features.cpp
    ${PROJECT_SOURCE_DIR}/src/audio/signal_process.cpp
    ${PROJECT_SOURCE_DIR}/src/image/color_map.cpp)
target_link_libraries(yuki
    Eigen3::Eigen
    fftw3
    ${FFMPEG_LIB}
    ${OpenCV_LIBS})
