cmake_minimum_required( VERSION 2.8)

project(ESIMATION)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib) 

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

include_directories(${PROJECT_SOURCE_DIR}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

link_directories(${PROJECT_SOURCE_DIR}/lib)

file(GLOB SRC_LSQ "./L*.cpp")

add_executable(lsq ${SRC_LSQ})

target_link_libraries(lsq vecmat gnss algorithm)

set(CMAKE_BUILD_TYPE "Debug")

add_executable(srif "./Srif.cpp")
target_link_libraries(srif vecmat gnss algorithm)