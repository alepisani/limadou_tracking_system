# filepath: c:\Users\user\Desktop\stageINFN\sim_mc\CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

project(sim_MC)

set(CMAKE_CXX_STANDARD 17)

# Handle path differences between Windows and WSL
if(CMAKE_HOST_UNIX AND EXISTS "/mnt/c/")
    # We're in WSL
    string(REGEX REPLACE "^/mnt/c/" "C:/" WIN_PATH "${CMAKE_SOURCE_DIR}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
    # We're in Windows
    set(WIN_PATH "${CMAKE_SOURCE_DIR}")
endif()

# Aggiungi i file sorgente
set(SOURCES
    src/LTrackerTrack.cpp
    src/stats.cpp
    src/display.cpp
    src/chip.cpp
    src/main.cpp
    src/eventdata.cpp
    src/LTrackerCluster.cpp
    src/simulations.cpp
)

# Aggiungi la directory include
include_directories(include)

# Linka le librerie ROOT
find_package(ROOT REQUIRED COMPONENTS Hist Graf Graf3d MathCore RIO)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Crea l'eseguibile
add_executable(sim_MC ${SOURCES})

# Linka le librerie ROOT
target_link_libraries(sim_MC ${ROOT_LIBRARIES})