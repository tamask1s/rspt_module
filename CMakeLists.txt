cmake_minimum_required(VERSION 3.4...3.18)
project(rspt_module)

set(CMAKE_CXX_STANDARD 14)

find_package(pybind11 REQUIRED)
find_package(Python COMPONENTS NumPy REQUIRED)

pybind11_add_module(rspt_module rspt_module.cpp)
target_link_libraries(rspt_module PRIVATE Python::NumPy)
