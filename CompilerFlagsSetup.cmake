# enable C++11
if (CMAKE_VERSION VERSION_GREATER 3.0.99)
    # let cmake set the flags
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
    # manually enable flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Enable compiler warnings and O2 optimization on gcc and g++
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse3 -march=native -fopenmp -Wall -Wno-unknown-pragmas -Wno-deprecated-declarations")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -msse3 -march=native -fopenmp -Wall -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -msse3 -march=native -Wall -Wno-unknown-pragmas -Wno-deprecated-declarations")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -msse3 -march=native -Wall -Wno-unknown-pragmas")
