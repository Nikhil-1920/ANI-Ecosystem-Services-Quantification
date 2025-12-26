cmake_minimum_required(VERSION 3.10)
project(EcoServiceANI)
set(CMAKE_CXX_STANDARD 17)

# --- MAC SPECIFIC ---
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/gdal;/opt/homebrew")
set(GDAL_INCLUDE_DIR "/opt/homebrew/include")
set(GDAL_LIBRARY "/opt/homebrew/lib/libgdal.dylib")

# Now find the package
find_package(GDAL REQUIRED)
include_directories(${GDAL_INCLUDE_DIR})

# Create the executable
add_executable(EcoServiceANI main.cpp)

# Link the library
target_link_libraries(EcoServiceANI ${GDAL_LIBRARY})