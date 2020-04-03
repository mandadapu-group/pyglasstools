set(SUBCOMPONENT_NAME potential)

set(_${SUBCOMPONENT_NAME}_sources
    module-potential.cpp
    )

pybind11_add_module(_${SUBCOMPONENT_NAME} SHARED ${_${SUBCOMPONENT_NAME}_sources} NO_EXTRAS)

# install the library
#install(TARGETS _${SUBCOMPONENT_NAME} 
#        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}/${SUBCOMPONENT_NAME}
#        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        PRE_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy pyglasstools/${SUBCOMPONENT_NAME}/${file}"
    )
ENDMACRO(copy_file)

set(files   __init__.py
    )

#install(FILES ${files}
#        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}/${SUBCOMPONENT_NAME}
#       )


foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_${SUBCOMPONENT_NAME} ALL DEPENDS ${files})
