CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#It is recommended that you copy all CMake related files in a local folder of yours,
#else whenever you pull/push there will be conflicts to resolve

SET(PROJECT_NAME "Cpp1")

project(${PROJECT_NAME})
# SET PATHS #
SET(PROJECT_DESTINATION "/home/dimihayl/Software/Python/Cpp1")
SET(GSL_INCLUDE "/home/dimihayl/Software/Python/Cpp1")
SET(GSL_LIB "/home/dimihayl/Software/Python/Cpp1")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DESTINATION}/bin)

add_executable(${PROJECT_NAME} ${PROJECT_DESTINATION}/main.cpp
)
# DO NOT CHANGE THE REST #

#string(APPEND CFLAGS " -O2 -std=c++11 -lgomp -pthread -fopenmp  ")
set(CFLAGS " -O2 -std=c++11 -pthread ${CFLAGS}")
#string(APPEND CFLAGS " -O2 -std=c++11 -pthread ")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${CFLAGS})

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})

target_link_libraries(${PROJECT_NAME} -L${GSL_LIB} -lgsl -lgslcblas)

