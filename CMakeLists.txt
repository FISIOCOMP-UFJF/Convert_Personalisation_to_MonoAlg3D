cmake_minimum_required(VERSION 2.8)

PROJECT(Converter)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fno-strict-aliasing -Wall -Wno-unused-function")

add_executable(Converter MACOSX_BUNDLE src/main.cpp )

target_link_libraries(Converter ${VTK_LIBRARIES})

add_subdirectory(src/calibration_library)

add_dependencies(Converter calibration)

target_link_libraries(Converter calibration dl m)


