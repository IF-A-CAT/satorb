cmake_minimum_required(VERSION 2.8)

add_library(vecmat SHARED "./VecMat.cpp")

add_library(algorithm SHARED "./Algorithm.cpp")

add_library(gnss SHARED "./GNSS.cpp")

target_link_libraries(gnss vecmat)

target_link_libraries(algorithm vecmat)



