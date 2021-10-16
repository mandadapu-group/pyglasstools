set(SUBCOMPONENT_NAME elasticitypetsc)

set(_${SUBCOMPONENT_NAME}_sources
    module-elasticity-petsc.cpp
    )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
find_package(MPI REQUIRED)

include_directories(${PETSC_INC} ${SLEPC_INC} ${MPI_CXX_INCLUDE_PATH})
pybind11_add_module(_${SUBCOMPONENT_NAME} SHARED ${_${SUBCOMPONENT_NAME}_sources} NO_EXTRAS)
target_link_libraries(_${SUBCOMPONENT_NAME} PRIVATE ${PETSC_LIB} ${SLEPC_LIB} ${MPI_CXX_LIBRARIES})
#set_target_properties(_${SUBCOMPONENT_NAME} PROPERTIES INSTALL_RPATH "\$ORIGIN/..;\$ORIGIN")

#install(TARGETS _${SUBCOMPONENT_NAME} 
#        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}/${SUBCOMPONENT_NAME}
#        )

set(files   __init__.py
    )

#install(FILES ${files}
#        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}/${SUBCOMPONENT_NAME}
#       )

foreach(file ${files})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()

add_custom_target(copy_${SUBCOMPONENT_NAME} ALL DEPENDS ${files})
