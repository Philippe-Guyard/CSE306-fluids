cmake_minimum_required(VERSION 3.5)

project(fluid_sim)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# target_include_directories(fluid_sim PRIVATE fluids)

set(SOURCES main.cpp
            fluids/vector.h
            fluids/vector.cpp
            fluids/poly.hpp
            fluids/diagram.hpp
            fluids/lbfgs.c
            fluids/lbfgs.h
            fluids/arithmetic_ansi.h
            fluids/fluid.hpp
            stb_image_write.h
            stb_image.h)
 
add_executable(fluid_sim ${SOURCES})