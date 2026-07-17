set(SOURCE_FILES tunecv.h tunecv.cpp)

add_library(calibration STATIC ${SOURCE_FILES})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(calibration ${VTK_LIBRARIES})

if(NOT MSVC)
    set_property(TARGET calibration PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()