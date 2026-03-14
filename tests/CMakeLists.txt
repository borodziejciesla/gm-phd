################################
# GTest
################################
project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        52eb810 # release-1.17.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

################################
# Tests
################################
add_executable(${TEST_PROJECT_NAME}
    calibrated_object_test.cpp
    gm_phd_test.cpp

    $<TARGET_OBJECTS:gm-phd>
)

target_link_libraries(${TEST_PROJECT_NAME}
    gtest 
    gtest_main
        
    gm-phd
)