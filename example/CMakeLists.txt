##########################################################
##################### GN-PHD Example #####################
##########################################################
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

add_executable(${EXAMPLE_PROJECT_NAME}
  cv_gm_phd.cpp
  example.cpp
  trajectory_generation.cpp
)

target_include_directories(${EXAMPLE_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:gm-phd,INTERFACE_INCLUDE_DIRECTORIES>
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS}
)

target_link_libraries(${EXAMPLE_PROJECT_NAME}
  gm-phd
  Python3::Python
)