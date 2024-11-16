cmake_minimum_required(VERSION 3.28)
project(MCalcTest)

set(CMAKE_CXX_STANDARD 14)

add_executable(MCalcTest main.cpp
        Matrix.cpp
        Matrix.h
)
