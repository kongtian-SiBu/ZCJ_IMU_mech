cmake_minimum_required(VERSION 3.10)

project(imu_homework)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

# Add head files for the target
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add source file for the target
set(SOURCE_FILE src/main.cpp)
add_executable(${PROJECT_NAME} ${SOURCE_FILE})

# Eigen
include_directories(/home/zhangnav/GNSSINS_CorCpp/common_libraries/eigen-3.4.0)