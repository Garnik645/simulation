cmake_minimum_required(VERSION 3.22.1)
project(stag)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -mtune=native -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native -ffast-math")

add_executable(stag main.c cJSON.c)
target_link_libraries("${PROJECT_NAME}" PRIVATE Threads::Threads)
