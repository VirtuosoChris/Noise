cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the project name and version
project(NoiseTests VERSION 1.0)

# Add the executable
add_executable(NoiseTests tests.cpp)

target_include_directories(NoiseTests PRIVATE ${CMAKE_SOURCE_DIR}/eigen)
