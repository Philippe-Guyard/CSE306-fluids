cmake_minimum_required(VERSION 3.5)

project(main)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCES main.cpp
            mesh_reader.hpp
            vector.cpp
            vector.h)

add_executable(main ${SOURCES})