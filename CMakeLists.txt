cmake_minimum_required(VERSION 2.9)

list(APPEND CMAKE_PREFIX_PATH "/usr/local/Qt-5.15.2")
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED)

project(FFT)

add_executable(fft fft.h fft.cpp test_fft.cpp)

add_executable(contrastePoisson fft.h fft.cpp
                                poisson.h poisson.cpp
                                contrastePoisson.cpp)
ImagineUseModules(contrastePoisson Images)

add_executable(clonePoisson fft.h fft.cpp
                            poisson.h poisson.cpp
                            clonePoisson.cpp)
ImagineUseModules(clonePoisson Images)
