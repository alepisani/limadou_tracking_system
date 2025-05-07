# filepath: c:\Users\user\Desktop\stageINFN\sim_mc\CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

project(sim_MC)

set(CMAKE_CXX_STANDARD 17)

# Aggiungi i file sorgente
set(SOURCES
    src/LTrackerTrack.cpp
    src/stats.cpp
    src/display.cpp
    src/main.cpp
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