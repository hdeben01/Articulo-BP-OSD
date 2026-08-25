# LINUX
cmake_minimum_required(VERSION 3.18)
project(OSD_decoder)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Buscar Python y pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Incluir Eigen
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Crear módulo
pybind11_add_module(OSD_decoder osd_pybind.cpp)

# Opciones de compilación Linux (opcional)
target_compile_options(OSD_decoder PRIVATE -O3 -Wall)



# WINDOWS
#cmake_minimum_required(VERSION 3.18)
#project(GaussianBP)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Buscar Python y pybind11
#find_package(Python COMPONENTS Interpreter Development REQUIRED)
#find_package(pybind11 REQUIRED)

# Incluir Eigen
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Crear módulo
#pybind11_add_module(GaussianBP GaussianBP.cpp)

# Opciones para UTF-8 (Python)
#add_compile_options(/utf-8)

