set(COMPONENT_NAME pyglasstools)

#### Add C++ Source Files
set(_${COMPONENT_NAME}_sources
    module-pyglasstools.cpp
    MPIInterface.cpp
    )

#### Finding External Packages

#Find EIGEN
list(APPEND INCLUDES extern/eigen/)

#Find Boost
find_package(Boost 1.50.0 REQUIRED serialization)
list(APPEND LIBRARIES ${Boost_LIBRARIES})
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

#Find GSL
find_package(GSL REQUIRED)   

#Find MPI
find_package(MPI REQUIRED)
list(APPEND LIBRARIES ${MPI_CXX_LIBRARIES})

#Find Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES extern/Aboria/src)
list(APPEND INCLUDES extern/Aboria/third-party)

#Find Spectra
list(APPEND INCLUDES extern/spectra/include)

#Cereal
list(APPEND INCLUDES extern/cereal/include)

#Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

#Add header files from the header-Only Libraries
include_directories(${INCLUDES} ${MPI_CXX_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

pybind11_add_module(_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} NO_EXTRAS)
target_link_libraries(_${COMPONENT_NAME} PRIVATE ${LIBRARIES})
#set_target_properties(_${COMPONENT_NAME} PROPERTIES INSTALL_RPATH "\$ORIGIN/..")
set_target_properties(_${COMPONENT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(_${COMPONENT_NAME} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

#install(TARGETS _${COMPONENT_NAME} 
#        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
#        )

add_custom_target(copy_${COMPONENT_NAME} ALL DEPENDS ${files})

#Link our .py scripts so that CMake updates them for newer builds
set(files   __init__.py
            utils.py
            CorrelationAnalysis.py
    )

#install(FILES ${files}
#        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
#       )

foreach(file ${files})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()

####Add our sub-modules
add_subdirectory(cgfunc)
add_subdirectory(potential)
add_subdirectory(observables)
add_subdirectory(nonaffine)
