# ==============================================================================
# Project Configuration
# ==============================================================================
cmake_minimum_required (VERSION 3.10)

# Set project name
project ("Ordered NR Tree")

set (CMAKE_CXX_STANDARD 17)

# ==============================================================================
# Include & Source
# ==============================================================================
include_directories ("${CMAKE_SOURCE_DIR}/include")
file(GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# ==============================================================================
# Build Target
# ==============================================================================
add_executable (main ${SOURCE_FILES})

# ======================================================================
# Runtime config copy
# ======================================================================
add_custom_target(copy_resources
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/config"
        "${CMAKE_CURRENT_BINARY_DIR}/config"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/data"
        "${CMAKE_CURRENT_BINARY_DIR}/data"
    COMMENT "Auto-copying resources to build folder..."
)

add_dependencies(main copy_resources)