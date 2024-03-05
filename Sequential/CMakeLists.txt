cmake_minimum_required(VERSION 3.10.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -fverbose-asm -masm=intel")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -c -Wa,-alh -masm=intel -fverbose-asm")

#=== Set name of our project
project(project2)

# Bring headers into the project
include_directories(include)
include_directories(include/H2)

# Add Sources
file(GLOB SOURCES "src/*")

set( LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/lib" )
add_library(project2 SHARED ${SOURCES})

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/scripts")
add_executable( main scripts/main.cxx )
target_link_libraries(main project2)

add_executable( vqmc scripts/optimized_vqmc.cxx )
target_link_libraries(vqmc project2)

add_executable( test scripts/benchmark.cxx )
target_link_libraries(test project2)
