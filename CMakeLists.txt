cmake_minimum_required(VERSION 3.1)
project(nhevo VERSION 1.0)

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
find_package(Eigen3 REQUIRED)
find_package(dlib REQUIRED)

set(HEADERS
  include/nhevo/solver.hh
  include/nhevo/solver_s3c2.hh
  include/nhevo/solver_s5c4.hh
  include/nhevo/solver_s7c6.hh
  )

set(SOURCES
  src/Sturm.cpp
  src/solver.cc
  src/solver_s3c2.cc
  src/solver_s5c4.cc
  src/solver_s7c6.cc
  )

add_library(solver SHARED
  ${HEADERS}
  ${SOURCES}
  )

target_include_directories(solver PUBLIC
  ${EIGEN3_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_link_libraries(solver PRIVATE
  dlib::dlib
  )

add_library(nhevo::solver ALIAS solver)


# -------------------
# install
# -------------------
include(GNUInstallDirs)

install(
  TARGETS solver
  EXPORT  nhevoTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nhevo
  )

install(
  EXPORT nhevoTargets
  FILE   nhevoTargets.cmake
  NAMESPACE nhevo::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nhevo
  )


# -------------------
# configure
# -------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/nhevoConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/nhevoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/nhevoConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nhevo
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/nhevoConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/nhevoConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nhevo
)

export(
  EXPORT nhevoTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/nhevoTargets.cmake
  NAMESPACE nhevo::
)

