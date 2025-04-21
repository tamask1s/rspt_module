cmake_minimum_required(VERSION 3.12)
project(rspt_module)

set(CMAKE_CXX_STANDARD 14)

find_package(pybind11 REQUIRED)

pybind11_add_module(rspt_module src/rspt_module.cpp)
