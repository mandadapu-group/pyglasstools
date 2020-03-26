set(COMPONENT_NAME pyglasstools)

set(_${COMPONENT_NAME}_sources
    module-pyglasstools.cpp
    )

#Find EIGEN
find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS})

#Find GSL
find_package(GSL REQUIRED)   

#Find Aboria
set(Aboria_LOG_LEVEL 2 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND ABORIA_INCLUDE_DIRS ./extern/Aboria/src)
list(APPEND ABORIA_INCLUDE_DIRS ./extern/Aboria/third-party)

include_directories(${EIGEN3_INCLUDE_DIRS} ${ABORIA_INCLUDE_DIRS})

pybind11_add_module(_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} NO_EXTRAS)

# install the library
install(TARGETS _${COMPONENT_NAME} 
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
        )
target_link_libraries(_${COMPONENT_NAME} GSL::gsl GSL::gslcblas)

    ################ Python only modules
# copy python modules to the build directory to make it a working python package
MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        PRE_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy pyglasstools/${COMPONENT_NAME}/${file}"
    )
ENDMACRO(copy_file)

set(files   __init__.py
            utils.py
            potential.py
            cgfunc.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_${COMPONENT_NAME} ALL DEPENDS ${files})
