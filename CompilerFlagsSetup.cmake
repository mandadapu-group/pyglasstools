# enable C++11
if (CMAKE_VERSION VERSION_GREATER 3.0.99)
    # let cmake set the flags
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
    # manually enable flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Enable compiler warnings and Ofast opeimization on gcc and g++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG -fopenmp -flto -march=native -Wall -Wno-unknown-pragmas -Wno-deprecated-declarations")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -DNDEBUG -fopenmp -flto -march=native -Wall -Wno-unknown-pragmas")
