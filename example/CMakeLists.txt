cmake_minimum_required(VERSION 3.1)
project(nhevo_ros_example VERSION 1.0)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()
message("Build type: " ${CMAKE_BUILD_TYPE})

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
  set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -Wall -Wextra -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# -------------------
# build
# -------------------
find_package(catkin REQUIRED COMPONENTS
  prophesee_event_msgs
  geometry_msgs
  roscpp
  )
find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(nhevo REQUIRED)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  )

add_executable(nhevo
  src/main.cpp
  src/nhevo.cc
  )

target_link_libraries(nhevo
  ${catkin_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  nhevo::solver
  )

catkin_package(
  INCLUDE_DIRS include
)
