cmake_minimum_required(VERSION 3.22)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 23)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

# Set the project
project(
  many-body-dynamics
  VERSION 0.0.0
  DESCRIPTION "Symbolic Calculations for Quantum Many-body Dynamics"
  LANGUAGES CXX)

# Add the source directory
add_subdirectory(src)
