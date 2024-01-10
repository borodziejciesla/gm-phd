# find_package(Python3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development)
if (${Python3_FOUND})
  include_directories(${Python3_INCLUDE_DIRS})
  message(VERBOSE ${Python3_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Python3 not found, please install it.")
endif()

##########################################################
##################### GM-PHD Example #####################
##########################################################
add_executable(${EXAMPLE_PROJECT_NAME}
  tg_example.cpp
  trajectory_generation.cpp
)

include_directories(${EXAMPLE_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:gm-phd,INTERFACE_INCLUDE_DIRECTORIES>
  ${Python3_INCLUDE_DIRS}
)

target_link_libraries(${EXAMPLE_PROJECT_NAME}
  gm-phd
  ${Python3_LIBRARIES}
)

##########################################################
################### GM-PHD-RHM Example ###################
##########################################################
add_executable(${EXAMPLE_RHM_PROJECT_NAME}
  rhm_example_matlab.cpp
  trajectory_generation.cpp
  sensor_data_reader.cpp
)

add_subdirectory(csv-parser)

include_directories(${EXAMPLE_RHM_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:gm-phd,INTERFACE_INCLUDE_DIRECTORIES>
  csv-parser/single_include
  ${Python3_INCLUDE_DIRS}
)

target_link_libraries(${EXAMPLE_RHM_PROJECT_NAME}
  gm-phd
  csv
  ${Python3_LIBRARIES}
)