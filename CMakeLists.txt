cmake_minimum_required(VERSION 3.24)
project(uav_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_compile_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
add_link_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)

set(UAV_WARNINGS -Wall -Wextra -Wpedantic)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(raylib REQUIRED)


add_library(uav_core STATIC
  src/frames.cpp
)
target_include_directories(uav_core PUBLIC include)
target_link_libraries(uav_core PUBLIC Eigen3::Eigen)
target_compile_options(uav_core PRIVATE ${UAV_WARNINGS})

add_library(uav_render STATIC
  src/render/render.cpp
)
target_include_directories(uav_render PUBLIC src)
target_link_libraries(uav_render PUBLIC uav_core raylib)
target_compile_options(uav_render PRIVATE ${UAV_WARNINGS})

add_executable(uav_sim src/main.cpp)
target_link_libraries(uav_sim PRIVATE uav_render)
target_compile_options(uav_sim PRIVATE ${UAV_WARNINGS})
