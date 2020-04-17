set(SUBCOMPONENT_NAME nonaffine)

set(_${SUBCOMPONENT_NAME}_sources
    module-nonaffine.cpp
    )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
find_package(MPI REQUIRED)

include_directories(${PETSC_INC} ${SLEPC_INC} ${MPI_CXX_INCLUDE_PATH})

pybind11_add_module(_${SUBCOMPONENT_NAME} SHARED ${_${SUBCOMPONENT_NAME}_sources} NO_EXTRAS)
target_link_libraries(_${SUBCOMPONENT_NAME} ${PETSC_LIB} ${SLEPC_LIB} ${MPI_CXX_LIBRARIES})

set(files   __init__.py
    )

foreach(file ${files})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()

add_custom_target(copy_${SUBCOMPONENT_NAME} ALL DEPENDS ${files})
