cmake_minimum_required(VERSION 3.14)
project(TCG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED COMPONENTS core imgproc imgcodecs ximgproc)

add_executable(TCG
  src/main.cpp
  src/load_edg.cpp
  src/load_cem.cpp
  src/contour_breaker.cpp
  src/image_gradient.cpp
  src/dp_gap.cpp
  src/contour_fill_gaps.cpp
  src/fac_graph.cpp
  src/break_t_junctions.cpp
  src/prune_noise.cpp
  src/merge_geom.cpp
  src/classify_bp.cpp
)

target_include_directories(TCG PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include" ${OpenCV_INCLUDE_DIRS})
target_link_libraries(TCG PRIVATE ${OpenCV_LIBS})

if(MSVC)
  target_compile_definitions(TCG PRIVATE _USE_MATH_DEFINES)
else()
  target_compile_definitions(TCG PRIVATE _DEFAULT_SOURCE)
endif()
