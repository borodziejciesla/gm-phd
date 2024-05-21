################################
# GTest
################################
project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        8a6feab # release-1.10.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

################################
# Tests
################################
add_executable(${TEST_PROJECT_NAME}
    hypothesis_test.cpp
)

target_include_directories(${TEST_PROJECT_NAME}
  PUBLIC
    gm-phd
)

target_link_libraries(${TEST_PROJECT_NAME}
    gtest 
    gtest_main

    Eigen3::Eigen
)