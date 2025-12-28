# ==============================================================================
# Project Configuration
# ==============================================================================
cmake_minimum_required (VERSION 3.8)

# Set project name
project ("Ordered NR Tree")

# Set C++17 standard (as specified in c_cpp_properties)
set (CMAKE_CXX_STANDARD 17)

# ==============================================================================
# Include Directories
# ==============================================================================
# Declare directory containing .h files (replaces -I in tasks.json)
include_directories ("${CMAKE_SOURCE_DIR}/include")

# ==============================================================================
# Source Files
# ==============================================================================
# Find all .cpp files (replaces *.cpp args in tasks.json)
file(GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# ==============================================================================
# Build Target
# ==============================================================================
# Create executable
add_executable (main ${SOURCE_FILES})

# ======================================================================
# Runtime config copy (IMPORTANT)
# ======================================================================
add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/config"
        "$<TARGET_FILE_DIR:main>/config"
    COMMENT "Copying config folder to output directory..."
)

# 2. Copy thư mục data
add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:main>/data"
    COMMENT "Copying data folder to output directory..."
)