project(pyglasstools)
cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include(SetupPybind11.cmake)
include(CompilerFlagsSetup.cmake)


# plugins must be built as shared libraries
if (ENABLE_STATIC)
    message(SEND_ERROR "Needs to be built as shared library!")
endif (ENABLE_STATIC)

set(BUILD_SHARED_LIBS on)

set(PYTHON_MODULE_BASE_DIR ${CMAKE_INSTALL_PREFIX})
message(STATUS "Install module to: " ${PYTHON_MODULE_BASE_DIR})

# add subdirectories
add_subdirectory(${PROJECT_NAME})
