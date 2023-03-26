cmake_minimum_required (VERSION 3.22)

project("gm-phd" VERSION 0.0.1
  DESCRIPTION "C++ implementation of GM-PHD Multiple Object Tracker"
  HOMEPAGE_URL "https://github.com/borodziejciesla/gm-phd")

set(TEST_PROJECT_NAME ${PROJECT_NAME}_tests)
set(EXAMPLE_PROJECT_NAME ${PROJECT_NAME}_example)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_EXAMPLE "Build example" OFF)

###############################################
################## Library ####################
###############################################
include(GNUInstallDirs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME}
  src/gm_phd_cv_pose.cpp
)

target_include_directories(${PROJECT_NAME}
  PUBLIC
    include
    ${EIGEN_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
  Eigen3::Eigen
)

###############################################
################# Examples ####################
###############################################
if(BUILD_EXAMPLE)
    project(${EXAMPLE_PROJECT_NAME})
    
    add_subdirectory(example)
endif() # BUILD_EXAMPLE

# ###############################################
# ################ Unit Tests ###################
# ###############################################
# if(BUILD_TESTS)    
#     project(${TEST_PROJECT_NAME})

#     add_subdirectory(tests)    
#     enable_testing()
#     add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})    
# endif() #   BUILD_TESTS