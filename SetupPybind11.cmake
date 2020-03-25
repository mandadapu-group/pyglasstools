# use pybind11's tools to find python and create python modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/pyglasstools/extern/pybind11/tools)

set(PYBIND11_PYTHON_VERSION 3)

# trick pybind11 tools to allow us to manage C++ flags
# cmake ignores this in 2.8, but when pybind11 sees this

#pybind11 has a cmake file to configure compatible c++ flags
include(pybind11Tools)

if (PYTHON_VERSION_MAJOR LESS 3)
    message(SEND_ERROR "PyGlassTools requires Python 3.x")
endif()
