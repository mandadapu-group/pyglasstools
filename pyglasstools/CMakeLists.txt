set(COMPONENT_NAME pyglasstools)

set(_${COMPONENT_NAME}_sources
    module-pyglasstools.cpp
    )

#Find EIGEN
find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS})
message( ${EIGEN3_INCLUDE_DIRS} )
               
pybind11_add_module(_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} NO_EXTRAS)

# install the library
install(TARGETS _${COMPONENT_NAME} 
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
        )

    ################ Python only modules
# copy python modules to the build directory to make it a working python package
MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy pyglasstools/${COMPONENT_NAME}/${file}"
    )
ENDMACRO(copy_file)

set(files   __init__.py
            utils.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_${COMPONENT_NAME} ALL DEPENDS ${files})