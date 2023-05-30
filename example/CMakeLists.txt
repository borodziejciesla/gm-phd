##########################################################
##################### GN-PHD Example #####################
##########################################################
find_package(Python3 COMPONENTS Interpreter Development)
if (${Python3_FOUND})
  include_directories(${Python3_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Python3 not found, please install it.")
endif()

add_executable(${EXAMPLE_PROJECT_NAME}
  tg_example.cpp
  trajectory_generation.cpp
)

include_directories(${EXAMPLE_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:gm-phd,INTERFACE_INCLUDE_DIRECTORIES>
#   $<TARGET_PROPERTY:eigen,INTERFACE_INCLUDE_DIRECTORIES>
  ${Python3_INCLUDE_DIRS}
)

target_link_libraries(${EXAMPLE_PROJECT_NAME}
  gm-phd
#   eigen
  ${Python3_LIBRARIES}
)