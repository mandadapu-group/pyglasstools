set(COMPONENT_NAME pyglasstools)

#### Add C++ Source Files
set(_${COMPONENT_NAME}_sources
    module-pyglasstools.cpp
    )

#### Finding External Packages

#Find EIGEN
#find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS})

#Find GSL
find_package(GSL REQUIRED)   

#Find Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND ABORIA_INCLUDE_DIRS extern/Aboria/src)
list(APPEND ABORIA_INCLUDE_DIRS extern/Aboria/third-party)

##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

#Add header files from the header-Only Libraries
include_directories(${EIGEN3_INCLUDE_DIRS} ${ABORIA_INCLUDE_DIRS} extern/spectra/include)

pybind11_add_module(_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} NO_EXTRAS)

#Link our .py scripts so that CMake updates them for newer builds
set(files   __init__.py
            utils.py
            calculator.py
    )
foreach(file ${files})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()

add_custom_target(copy_${COMPONENT_NAME} ALL DEPENDS ${files})

####Add our sub-modules
add_subdirectory(cgfunc)
add_subdirectory(potential)
add_subdirectory(observables)
add_subdirectory(nonaffine)
