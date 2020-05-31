set(SUBCOMPONENT_NAME cgfunc)

#### Add C++ Source Files
set(_${SUBCOMPONENT_NAME}_sources
    module-cgfunc.cpp
    )

pybind11_add_module(_${SUBCOMPONENT_NAME} SHARED ${_${SUBCOMPONENT_NAME}_sources} NO_EXTRAS)

#### Link the module with GSL libraries
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(_${SUBCOMPONENT_NAME} GSL::gsl GSL::gslcblas)

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
