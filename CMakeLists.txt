cmake_minimum_required (VERSION 3.22)

project("gm-phd" VERSION 0.0.1
  DESCRIPTION "C++ implementation of GM-PHD Multiple Object Tracker"
  HOMEPAGE_URL "https://github.com/borodziejciesla/gm-phd")

set(TEST_PROJECT_NAME ${PROJECT_NAME}_tests)
set(EXAMPLE_PROJECT_NAME ${PROJECT_NAME}_example)

option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLE "Build example" OFF)

###############################################
################## Library ####################
###############################################
include(GNUInstallDirs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(TBB REQUIRED)

add_library(${PROJECT_NAME}
  src/cv_motion_model.cpp
)

if(CMAKE_VERSION VERSION_LESS "3.25")
    target_compile_options(${PROJECT_NAME} PUBLIC -std=c++2c)
else()
    target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_26)
endif()

target_include_directories(${PROJECT_NAME}
  PUBLIC
    include
    ${EIGEN_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
  Eigen3::Eigen
  TBB::tbb
)

###############################################
################# Examples ####################
###############################################
if(BUILD_EXAMPLE)    
    add_subdirectory(example)
endif() # BUILD_EXAMPLE

# ###############################################
# ################ Unit Tests ###################
# ###############################################
if(BUILD_TESTS)    
    add_subdirectory(tests)    
    enable_testing()
    add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})    
endif() #   BUILD_TESTS